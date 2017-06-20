#!/usr/bin/env python
#python 2.7.5 requires biopython
#TE-Reunite.py
#Version 0.0.1 Adam Taranto, June 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

###################################################################################
# For a set of reference transposons, collect instances from one or more genomes. #
# Use aligned collections to build HMMs capturing TE family diversity.            #
###################################################################################

import sys
import os
import time
import argparse
from collections import Counter
from collections import namedtuple
from copy import deepcopy
from Bio import SeqIO
import Bio.Align.Applications
#import numpy as np
#import seaborn as sns
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
	#else: #Report label pairs
		#print("Setting Genome labels:")
		#for i in range(len(args.genomes)):
			#print(str(args.genomeLabels[i]) + "==" + args.genomes[i])

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
	checkAligners(args)
	# Make outDir if does not exist else set to current dir.
	outDir = checkPaths(args.outDir)
	if args.doAlign:
		# Make alnDir if does not exist else set to current dir.
		alnDir = checkPaths(args.alnDir)
	else:
		alnDir = None
	# Set temp directory
	tempDir = makeTemp()
	# Return paths
	return outDir,tempDir,alnDir

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

def writeClusters(allHits,refMaster,genMaster,outDir,tempDir,getAnchors=False):
	clusterOutPaths = list()
	# Open Error Log file
	# refMaster {refID:SeqRecord}
	# genMaster = {GenLabel:{ChromName:SeqRecord}}
	# allHits = {repeatName:[Hit1,Hit2,Hit3]}
	# ['refTE','ChrName','GenLabel','Strand','Gaps','Idt','HitLen','Qstart','Qend','HitStart','HitEnd']
	for refTE in refMaster.keys():
		outPath = os.path.join(outDir, "Cluster_" + refTE + ".fa")
		clusterOutPaths.append(outPath)
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
	return clusterOutPaths

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

"""
def buildAnchors(target,blastHits):
	anchors = list()
	anchorItem = namedtuple('Seq1','Seq2','Start1','Start2','hitlen','Score','ID')
	for hit in blastHits[target]:
		if hit.gaps == 0:
			# Set Name Target_GenID_ChrmID_Start_End_Strand
			anchors.append(anchorItem(target,hitName,Qstart,Hstart,hitLen,0,ID))
	# Sort anchors by ID
	# Update score for anchors (high ID to low)
	return anchors

def writeAnchors(tempPath,anchors):
	# Unpack and write to tempPath
	pass

def alignSet(outPath,tempPath):
	# If anchors not none
	# Write anchors to temp
	# Align with anchors
	# Else - Align without anchors
	pass

def countIntersects(clusters):
	# Get all unique target pairs
	# Count number of hits in target1 which overlap target2 hits
	# Plot heatmap with seaborn
	pass

# Optional: Attempt alignment of reunited repeat families
def alignSets():
	pass
# Optional: Produce heatmap of overlapping repeat hits  
def plotIntersects():
	pass
"""

def main(args):
	# Check for required files + make output directories
	outDir,tempDir,alnDir = dochecks(args)
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
	# Extract hits for each refRepeat, write clusters and return summary object 
	# list of tuples: [(clustPath,anchorPath,clustCount,anchorCount)]
	clusterPaths = writeClusters(allHits,refMaster,genMaster,outDir,tempDir,getAnchors=args.anchors)

	# Optional: Attempt alignment of reunited repeat families
	#if args.doAlign:
	#	alignSets(clusters=clusterPaths,outPath=alnDir,alnTool=args.AlnTool,outFormat=args.outAlnFormat)
	# Optional: Produce heatmap of overlapping repeat hits  
	#if args.heatmap:
	#	plotIntersects(allHits,refMaster,outDir)
	# Clean temp directory
	 # Delete tempDir

if __name__== '__main__':
	__version__ = '0.0.1'
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
	parser.add_argument("--heatmap",
								action="store_true",
								help="Plot heatmap for number of overlapping hits between all \
								pairs of reference repeats. Indicates potential duplicated or nested \
								reference sequences.")
	parser.add_argument("--anchors",
								action="store_true",
								help="Calculate alignment anchor table for all ungapped hits. \
								For use with DIALIGN.")

	args = parser.parse_args()
	main(args)