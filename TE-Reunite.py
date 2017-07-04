#!/usr/bin/env python
#python 2.7.5 requires biopython
#TE-Reunite.py
#Version 0.1.1 Adam Taranto, June 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

###################################################################################
# For a set of reference transposons, collect instances from one or more genomes. #
# Use aligned collections to build models capturing TE family diversity or        #
# produce deRIP'd reference sequence.                                             #
###################################################################################

import sys
import os
import time
import argparse
import itertools
from collections import Counter
from collections import namedtuple
from copy import deepcopy
from Bio import SeqIO
import Bio.Align.Applications
from _version import __version__

def makeTemp():
	return os.path.join(os.getcwd(),"_".join(str(time.time()).split(".") + ["temp"]))

def checkPaths(pathIs):
	if pathIs:
		fullpath = os.path.abspath(pathIs)
		if not os.path.isdir(fullpath):
			os.makedirs(fullpath)
	else:
		fullpath = os.getcwd()
	return fullpath

def checkUnique(somelist):
	IDcounts = Counter(somelist)
	duplicates = [k for k, v in IDcounts.items() if v > 1]
	if duplicates:
		print("Duplicated elements in arguements. Quiting.")
		print(duplicates)
		sys.exit(1)
	else:
		pass

def checkAligners(args):
	if args.doAlign:
		toolDict = {"MUSCLE":"alnMuscleCommandline",
					"clustalW":"ClustalwCommandline",
					"DIALIGN":"DialignCommandline"}
		if toolDict[args.AlnTool] not in dir(Bio.Align.Applications):
			raise IOError("%s is not available to biopython." % args.AlnTool)
		# Check DIALIGN is available
		if args.anchors and "DialignCommandline" not in dir(Bio.Align.Applications):
			print("Dialign not available to Biopython, alignment anchors with be generated but not used.")

def checkSameLen(args):
	# Check same number of genome labels and files
	if len(args.genomes) != len(args.genomeLabels):
		raise IOError('Different number of Labels and Genome files.')

def checkExists(args):
	# Check all genome files exist
	for gen in args.genomes:
		if not os.path.isfile(gen):
			raise IOError('File %s does not exist.' % gen)
	# Check BLAST TAB file exists if path provided
	if args.blastTAB and not os.path.isfile(args.blastTAB):
			raise IOError('File %s does not exist.' % args.blastTAB)
	# Check RM-OUT file exists if path provided
	if args.rmOut and not os.path.isfile(args.rmOut):
			raise IOError('File %s does not exist.' % args.rmOut)
	# Check repeat library exists
	if not os.path.isfile(args.repeats):
			raise IOError('File %s does not exist.' % args.repeats)

def dochecks(args):
	# Check files exist
	checkExists(args)
	# Check that genome files and labels are unique
	checkUnique(args.genomes)
	checkUnique(args.genomeLabels)
	# Check same number of genome labels and files
	checkSameLen(args)
	# Check aligner is installed
	#checkAligners(args)
	# Make outDir if does not exist else set to current dir.
	outDir = checkPaths(args.outDir)
	#if args.doAlign:
		# Make alnDir if does not exist else set to current dir.
	#	alnDir = checkPaths(args.alnDir)
	#else:
	#	alnDir = None
	# Set temp directory
	#tempDir = makeTemp()
	# Return paths
	return outDir#,tempDir#,alnDir

def chunkstring(string, length=80):
	return (string[0+i:length+i] for i in range(0, len(string), length))

def revComplement(seq):
	revcompl = lambda x: ''.join([ {'A':'T', 'T':'A', 'U':'A', 
									'G':'C', 'C':'G', 'Y':'R', 
									'R':'Y', 'S':'S', 'W':'W', 
									'K':'M', 'M':'K', 'B':'V', 
									'D':'H', 'H':'D', 'V':'B', 
									'N':'N'}[Z] for Z in x][::-1])
	return revcompl(seq)

def getStrand(start,end):
	if start <= end:
		return ("+",start,end)
	elif start > end:
		return ("C",end,start)

def importRefseqs(targetPath,stripHash=False):
	"""Populate dictionary with reference TEs, checking that names are unique."""
	refMaster = dict()
	if not stripHash:
		for rec in SeqIO.parse(targetPath, "fasta"):
			if rec.id not in refMaster.keys():
				refMaster[rec.id] = rec
			else:
				# If ref TE name occurs > 1, die.
				print("WARNING: Sequence name %s occurs multiple times in reference repeat fasta." % str(rec.id))
				sys.exit(1)
	else:
		for rec in SeqIO.parse(targetPath, "fasta"):
			if str(rec.id).split('#',1)[0] not in refMaster.keys():
				refMaster[str(rec.id).split('#',1)[0]] = rec
			else:
				# If ref TE name occurs > 1, die.
				print("WARNING: Sequence name %s occurs multiple times in reference repeat fasta." % str(rec.id).split('#',1)[0])
				sys.exit(1)
	return refMaster

def importGenomes(genomeList,GenLabels,oldSeq=None):
	"""Populate dictionary with master set of seq records."""
	# Init sequence object as nested dict = GenLabel:{ChromName:SeqRecord}
	if not oldSeq:
		SeqMaster = dict()
		for label in GenLabels:
			SeqMaster[label] = dict()
	else: 
		# Copy existing sequence dict and add new labels
		SeqMaster = deepcopy(oldSeq)
		for label in GenLabels:
			if label not in SeqMaster.keys():
				SeqMaster[label] = dict()
	# Init blank chrom to genome map
	mapC2G = dict()
	# Import sequences from each genome file
	for inFasta,label in zip(genomeList,GenLabels):
		for seq_record in SeqIO.parse(inFasta, "fasta"):
			SeqMaster[label][seq_record.id] = seq_record
			# If first instance of chromosome name, store with genome label
			if seq_record.id not in mapC2G.keys():
				mapC2G[seq_record.id] = label
			else:
				# If Chrom name occurs > 1 time across genomes, die.
				print("WARNING: Sequence name %s occurs multiple times." % str(seq_record.id))
				sys.exit(1)
	return SeqMaster,mapC2G

def importBLAST(infile=None,QuerySeq=None,chr2Gen=None,minCov=0.9,minID=0.9,eVal=0.001):
	"""Import vaild blast hits to dict with structure {repeatName:[Hit1,Hit2,Hit3]} 
	where Hits are named tuples containing relevant hit information."""
	""" BLASTn output format 6
	0		qseqid			query (e.g., gene) sequence id
	1		sseqid			subject (e.g., reference genome) sequence id
	2		pident			percentage of identical matches
	3		length			alignment length
	4		mismatch		number of mismatches
	5		gapopen			number of gap openings
	6		qstart			start of alignment in query
	7		qend			end of alignment in query
	8		sstart			start of alignment in subject
	9		send			end of alignment in subject
	10		evalue			expect value
	11		bitscore		bit score
	"""
	# Init dict
	blastHits = dict()
	# Add refTE names as keys
	for refTE in QuerySeq.keys():
		blastHits[refTE] = list()
	# Init named tuple structure
	hitItem = namedtuple("HitRecord",['refTE','ChrName','GenLabel',
										'Strand','Gaps','Idt','HitLen',
										'Qstart','Qend','HitStart',
										'HitEnd'])
	# Read in BLAST tab file
	with open(infile) as src:
		data = [rec.strip().split('\t') for rec in src]
		for line in data:
			# Ignore lines begining with '#'
			if line[0][0] == "#":
				continue
			#Check that hit chrom exists in imported genomes
			if str(line[1]) not in chr2Gen.keys():
				continue
			# Calculate proportion of query covered by hit alignment
			hitCoverage = int(line[3])/len(QuerySeq[str(line[0])].seq)
			# Ignore if eVal is above threshold, or ID or Coverage are too low 
			if float(line[10]) >= eVal or float(line[2])/100 <= minID or hitCoverage <= minCov:
				continue
			# Check strand and add hitItem to queryTE hitlist
			strand,hitS,hitE = getStrand(int(line[8]),int(line[9]))
			blastHits[str(line[0])].append(hitItem(refTE=line[0], ChrName=str(line[1]),
													GenLabel=chr2Gen[str(line[1])], Strand=strand,
													Gaps=int(line[5]), Idt=float(line[2])/100, HitLen=int(line[3]), 
													Qstart=int(line[6]), Qend=int(line[7]),
													HitStart=hitS, HitEnd=hitE
													))
			# Sort entries by Genome, Chrm, Start, End
			blastHits[str(line[0])] = sorted(blastHits[str(line[0])], key = lambda x: (x.GenLabel, x.ChrName, x.HitStart, x.HitEnd))
	return blastHits

def importRM(infile=None,QuerySeq=None,chr2Gen=None,minCov=0.9,minID=0.9,score=100,maxInsert=0.1,maxFlanks=0.1):
	"""Import vaild RM.out hits to dict with structure {repeatName:[Hit1,Hit2,Hit3]} 
	where Hits are named tuples containing relevant hit information."""
	""" RM out
	0		SW score
	1		% substitutions in matching region compared to the consensus (Divergence) [%SNPs]
	2		% of gaps in hit (deleted bp relative to repeat) [Calc = (#Gaps in hit) / alignmentLength] 
	3		% of gaps in reference repeat (inserted bp in hit, relative to reference repeat) [Calc = (#Gaps in repeat) / alignmentLength]
	4		name of sequence containing hit
	5		hitStart
	6		hitEnd
	7		no. of bases in query sequence past the ending position of match
	8		strand
	9		Name of repeat (query sequence)
	10		Class of the repeat (string in lib name following '#')
	11		no. of bases in (complement of) the repeat consensus sequence prior to beginning of the match (so 0 means that the match extended all the way to the end of the repeat consensus sequence)
	12		Starting position of match in database sequence (using top-strand numbering)
	13		Ending position of match in database sequence
	14		Sequential hit identifier	
	"""
	# Init dict
	rmHits = dict()
	# Add refTE names as keys
	for refTE in QuerySeq.keys():
		rmHits[refTE] = list()
	# Init named tuple structure
	hitItem = namedtuple("HitRecord",['refTE','ChrName','GenLabel',
										'Strand','Gaps','Idt','alignLen',
										'Qstart','Qend','HitStart',
										'HitEnd'])
	# Read in RM out file
	with open(infile) as src:
		data = [rec.strip().split() for rec in src]
		for line in data:
			if line: 
				if line[0][0].isdigit():
					# Calculate proportion of query covered by hit alignment
					if str(line[8]) == 'C':
						Qstart = int(line[13])
						Qend = int(line[12])
					else:
						Qstart = int(line[11])
						Qend = int(line[12])

					refLen = Qend - Qstart # Length of refTE included in alignment
					rawHitLen = int(line[6])-int(line[5]) # Includes any internal insertions, may be > than refLen
					hitCov = refLen/len(QuerySeq[str(line[9])].seq) # Proportion of refTE covered by alinged segments of hit, does not include insertions
					hitID = (100-float(line[1]))/100 # Prop. matches in aligned positions
					hitInsert = rawHitLen/refLen #Raw hit interval as proportion of alignment to reference (i.e. (Aligned regions of hit + insertions spaned by RM) / Aligned region of refTE, value > 1 indicates insertions)
					alignLen = int(round(refLen-1) * (1+(float(line[2])/100))) # Aligned region of refTE + gaps in hit
					gapsInHit = int(round(alignLen - (refLen))) # ~Positions deleted in hit relative to refTE
					gapsInRef = int(round(alignLen * (float(line[3])/100))) # ~Positions inserted in hit relative to refTE
					# Ignore if SW-score, ID or Coverage are too low 
					if float(line[0]) < score or hitID < minID or hitCov < minCov or hitInsert > (maxInsert+1) or hitCov > (maxFlanks+1):
						continue
					else:
						rmHits[str(line[9])].append(hitItem(refTE=str(line[9]), ChrName=str(line[4]),
																GenLabel=chr2Gen[str(line[4])], Strand=str(line[8]),
																Gaps=gapsInRef, Idt=hitID, alignLen=alignLen, 
																Qstart=Qstart, Qend=Qend,
																HitStart=int(line[5]), HitEnd=int(line[6])
																))
						rmHits[str(line[9])] = sorted(rmHits[str(line[9])], key = lambda x: (x.GenLabel, x.ChrName, x.HitStart, x.HitEnd))
	return rmHits

def writeClusters(allHits,refMaster,genMaster,outDir,WriteSummary=True,SkipZeros=False,writeSeqs=True):
	##clusterOutPaths = list()
	# Open Error Log file
	# refMaster {refID:SeqRecord}
	# genMaster = {GenLabel:{ChromName:SeqRecord}}
	# allHits = {repeatName:[Hit1,Hit2,Hit3]}
	# ['refTE','ChrName','GenLabel','Strand','Gaps','Idt','HitLen','Qstart','Qend','HitStart','HitEnd']
	for refTE in refMaster.keys():
		if len(allHits[refTE]) == 0 and SkipZeros:
			continue
		else:
			outPath = os.path.join(outDir, refTE + "_ReuniteHits_" + str(len(allHits[refTE])) +".fa")
			##clusterOutPaths.append(outPath)
			outHandle = open(outPath,'w')
			# Write reference TE
			outHandle.write(">%s \n" % refTE)
			for line in chunkstring(str(refMaster[refTE].seq)):
				outHandle.write(line + "\n")
			# Write hits to cluster
			for label,seq in getFragments(genMaster,allHits,refTE):
				outHandle.write(">%s \n" % label)
				for line in chunkstring(seq):
					outHandle.write(line + "\n")
			outHandle.close()
	##return clusterOutPaths

def getFragments(genMaster,allHits,refTE):
	for x in allHits[refTE]:
		if x.Strand == '+':
			txtStrand = 'F'
		else:
			txtStrand = 'C'
		#Compose fragment label
		lable = "_".join([x.refTE, x.ChrName, x.GenLabel, txtStrand, str(x.HitStart), str(x.HitEnd)])
		# Extract fragment, reverse complement if hit is on rev strand
		if x.Strand == "C":
			seq = revComplement(str(genMaster[x.GenLabel][x.ChrName].seq[x.HitStart:x.HitEnd]))
		else:
			seq = str(genMaster[x.GenLabel][x.ChrName].seq[x.HitStart:x.HitEnd])
		yield (lable,seq)

def getOverlap(a, b):
	return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def checkOverlap(Astart,Aend,Bstart,Bend,pCov=None):
	if (Aend >= Bstart and Aend <= Bend) or (Astart <= Bend and Astart >= Bstart) or (Astart >= Bstart and Aend <= Bend) or (Bstart >= Astart and Bend <= Aend):
		if not pCov:
			return True
		else:
			Alen = Aend - Astart
			Blen = Bend - Bstart
			insLen = getOverlap((Astart,Aend),(Bstart,Bend))
			if insLen/Alen >= pCov or insLen/Blen >= pCov:
				return True
			else:
				return False
	else:
		return False

def countIntersects(allHits,pCov):
	countDict = dict()
	for refA,refB in itertools.combinations(allHits.keys(), 2):
		for hitA in allHits[refA]:
			for hitB in allHits[refB]:
				if hitB.GenLabel == hitA.GenLabel and hitA.ChrName == hitB.ChrName:
					if checkOverlap(hitA.HitStart,hitA.HitEnd,hitB.HitStart,hitB.HitEnd,pCov):
						if frozenset([refA,refB]) not in countDict.keys():
							countDict[frozenset([refA,refB])] = 1
						else:
							countDict[frozenset([refA,refB])] += 1
	return countDict

def groupOverlaps(d):
	""" For all unique pairs of refTEs sharing one or more overlapping (of >=1bp) 
		hits in any target genome, merge pairs which share at least one refTE in 
		common. Report merged overlap clusters."""
	d = dict((k, v) for k, v in d.items() if v >= 1)
	outlist = list()
	inlist = sorted([sorted(list(x)) for x in d.keys()])
	if not inlist:
		pass
	else:
		outlist.append(inlist[0]) # Add the first item from input list to outlists
		for l in inlist[1:]: # Starting from second item in input list
			merge = False # Reset merge log
			listSet = set(l) # Convert to set
			for index in range(len(outlist)): # For each outlist
				rset = set(outlist[index]) # Convert to set
				if len(listSet & rset) != 0 and not merge: # If >=1 value shared value between inlist and current outlist
					outlist[index] = sorted(list(listSet | rset)) # Merge lists and update position in outlists
					merge = True # Log merge event
					mergeIdx = index
				elif len(listSet & rset) != 0 and merge:
					outlist[mergeIdx] = sorted(list(set(outlist[mergeIdx]) | rset))
					outlist[index] = []
			if not merge: # Add new list to outlist if unmerged
				outlist.append(l)
			outlist = [x for x in outlist if x != []]
	# Return list of clusters sorted from largest to smallest.
	return sorted(outlist, key=len,reverse=True)

def	writeOverlaps(outDir,clusters):
	outPath = os.path.join(outDir, "Hit_Overlap_Clusters" + ".txt")
	outHandle = open(outPath,'w')
	for i in range(len(clusters)):
		outHandle.write("Group_%s:\t" % str(i) + '\t'.join(clusters[i]) + "\n")
	outHandle.close()

def	writeRedundant(outDir,clusters,refMaster):
	"""clusters is a list of lists"""
	for i in range(len(clusters)): # For each group
		outPath = os.path.join(outDir, "Group_%s.fa" % str(i))
		outHandle = open(outPath,'w')
		for refTE in clusters[i]: # For each sequence in group
			outHandle.write(">%s \n" % refTE) # Write sequence name
			for line in chunkstring(str(refMaster[refTE].seq)):
				outHandle.write(line + "\n")
		outHandle.close()

def writeGenomeSummary(outDir, allHits, genMaster):
	tracker = dict()
	for refTE in allHits.keys():
		tracker[refTE] = dict()
		for gen in genMaster.keys():
			tracker[refTE][gen] = 0
		for hit in allHits[refTE]:
			tracker[refTE][hit.GenLabel] += 1

	GenNames = sorted(genMaster.keys())
	outPath = os.path.join(outDir, "Summary_Hits_by_Genome.txt")
	outHandle = open(outPath,'w')
	#Write header
	outHandle.write("Name\t" + '\t'.join(GenNames) + '\n')
	for x in sorted(tracker.keys()):
		countList = list()
		for name in GenNames:
			countList.append(str(tracker[x][name]))
		outHandle.write(str(x) + "\t" + '\t'.join(countList) + '\n')
	outHandle.close()

def main(args):
	# Check for required files + make output directories
	outDir = dochecks(args)
	# Import target genomes as {GenLabel:{ChromName:SeqRecord}}
	# Also return seq to genome map as {SeqName:GenLabel}
	genMaster,chr2Gen = importGenomes(args.genomes,args.genomeLabels)
	# Import reference repeat library, check for unique names.
	if args.blastTAB:
		refMaster = importRefseqs(args.repeats)
	elif args.rmOut:
		refMaster =importRefseqs(args.repeats,stripHash=True)
	# Import BLAST or RM hits, filter for query coverage and identity, determine orientation
	# Dict key by refName
	if args.blastTAB:
		allHits = importBLAST(infile=args.blastTAB, QuerySeq=refMaster,
								chr2Gen=chr2Gen, minCov=args.mincov,
								minID=args.minID, eVal=args.eVal)
	elif args.rmOut:
		allHits = importRM(infile=args.rmOut, QuerySeq=refMaster,
								chr2Gen=chr2Gen, minCov=args.mincov,
								minID=args.minID, score=args.minSW)
	else:
		sys.exit(1)

	# Find number of overlapping hits between any two refTEs
	if args.reportoverlaps:
		overlapDict = countIntersects(allHits,args.pOverlap)
		refClusters = groupOverlaps(overlapDict)
		writeOverlaps(outDir,refClusters)
		if args.redundantRefs:
			writeRedundant(outDir,refClusters,refMaster)

	# Extract hits for each refRepeat, write clusters
	writeClusters(allHits,refMaster,genMaster,outDir,SkipZeros=args.onlyhits)
	# Write summary of hits per refRepeat found per target genome
	writeGenomeSummary(outDir, allHits, genMaster)


if __name__== '__main__':
	__version__ = '0.1.1'
	###Argument handling.
	parser = argparse.ArgumentParser(
								description='For a set of reference transposons, collect instances from one or more genomes.',
								prog='TE-Reunite')
	parser.add_argument('-v', 
								'--version', 
								action='version', 
								version='%(prog)s {version}'.format(version=__version__))
	parser.add_argument("-r", 
								"--repeats",
								type=str,
								required=True,
								help="Fasta formated library of reference repeats. Note: Names must be unique."
								)
	parser.add_argument("-g", 
								"--genomes",
								nargs='+',
								required=True,
								help="Space delimited list of reference genomes used as BLAST targets. Fasta format. \
								Note: Names of sequences must be unique across ALL genomes provided."
								)
	parser.add_argument("-l", 
								"--genomeLabels",
								nargs='+',
								required=True,
								help="Space delimited list of genome labels. \
								Note: Lables must be unique and same length as genome files list."
								)
	parser.add_argument("-b", 
								"--blastTAB",
								type=str,
								default=None,
								help="BLAST output format 6."
								)
	parser.add_argument("--rmOut",
								type=str,
								default=None,
								help="Repeatmasker outfile."
								)
	parser.add_argument("-s", 
								"--minSW",
								type=float,
								default=0.9,
								help="Minimum Smith-Waterman alignment score if providing hits in RepeatMasker.out format."
								)
	parser.add_argument("-c", 
								"--mincov",
								type=float,
								default=0.9,
								help="Minimum coverage of reference sequence required for vaild hit. Range: 0-1"
								)
	parser.add_argument("-e", 
								"--eVal",
								type=float,
								default=0.001,
								help="Minimum e-Value required for vaild hit."
								)
	parser.add_argument("-i", 
								"--minID",
								type=float,
								default=0.9,
								help="Minimum sequence identity for vaild hit. Range: 0-1"
								)
	parser.add_argument("-o",
								"--outDir",
								type=str,
								default=None, 
								help="Write reunited transposon families to this directory.")
	"""
	parser.add_argument("--doAlign",
								action="store_true",
								help="If set attempt to align reunited repeat families.")
	parser.add_argument("--alnDir",
								type=str,
								default=None, 
								help="Write aligned transposon families to this directory.")
	parser.add_argument('--AlnTool',
								default="clustalW",
								choices=["clustalW","DIALIGN","MUSCLE"],
								help='Attempt to align clusters with this tool. \
								Note: Will first attempt alignment with DIALIGN if anchors set.')
	parser.add_argument('--outAlnFormat',
								default="fasta",
								choices=["clustal","emboss","fasta","fasta-m10","ig","nexus","phylip","phylip-sequential","phylip-relaxed","stockholm"],
								help='Optional: Write alignment including reference sequence to file of format X.')
	"""
	parser.add_argument("--onlyhits",
								action="store_true",
								help="If set, suppress output clusters containing only the reference repeat with no additional hits.")
	parser.add_argument("--pOverlap",
								type=float,
								default=0.8,
								help="Minimum overlap between hits from two refTEs for hit loci to be considered as shared. \
								Length of intersect as prop. of either hit in pair of hits must be >= to this value. \
								i.e. A small insert in a larger repeat will have an overlap of 1. Range: 0-1"
								)
	parser.add_argument("--reportoverlaps",
								action="store_true",
								help="Report clusters of reference repeats which share hit locations overlapping >= 1bp with at least \
								one other member of the cluster. Use as guide to curate and merge redundant reference repeats.")
	parser.add_argument("--redundantRefs",
								action="store_true",
								help="If set with 'reportoverlaps', print groups of potentially redundant reference repeats \
								to multi FASTA files for inspection.")


	args = parser.parse_args()
	main(args)