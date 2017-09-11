#!/usr/bin/env python
#python 2.7.5 requires biopython
#Adam Taranto, June 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

###################################################################################
# For a set of reference transposons, collect instances from one or more genomes. #
# Use aligned collections to build models capturing TE family diversity or        #
# produce deRIP'd reference sequence.                                             #
###################################################################################

import argparse
import reunite
import sys

def getArgs():
	###Argument handling.
	parser = argparse.ArgumentParser(description='For a set of reference transposons, collect instances from one or more genomes.', prog='reunite')
	parser.add_argument("-r",  "--repeats", type=str, required=True, help="Fasta formated library of reference repeats. Note: Names must be unique." )
	parser.add_argument("-g",  "--genomes", nargs='+', required=True, help="Space delimited list of reference genomes used as BLAST targets. Fasta format. \ Note: Names of sequences must be unique across ALL genomes provided." )
	parser.add_argument("-l",  "--genomeLabels", nargs='+', required=True, help="Space delimited list of genome labels. \ Note: Lables must be unique and same length as genome files list." )
	parser.add_argument("-b",  "--blastTAB", type=str, default=None, help="BLAST output format 6." )
	parser.add_argument("--rmOut", type=str, default=None, help="Repeatmasker outfile." )
	parser.add_argument("-s",  "--minSW", type=float, default=0.9, help="Minimum Smith-Waterman alignment score if providing hits in RepeatMasker.out format." )
	parser.add_argument("-c",  "--mincov", type=float, default=0.9, help="Minimum coverage of reference sequence required for vaild hit. Range: 0-1" )
	parser.add_argument("-e",  "--eVal", type=float, default=0.001, help="Minimum e-Value required for vaild hit." )
	parser.add_argument("-i",  "--minID", type=float, default=0.9, help="Minimum sequence identity for vaild hit. Range: 0-1" )
	parser.add_argument("-o", "--outDir", type=str, default=None,  help="Write reunited transposon families to this directory.")
	parser.add_argument("--onlyhits", action="store_true", help="If set, suppress output clusters containing only the reference repeat with no additional hits.")
	parser.add_argument("--pOverlap", type=float, default=0.8, help="Minimum overlap between hits from two refTEs for hit loci to be considered as shared. \ Length of intersect as prop. of either hit in pair of hits must be >= to this value. \ i.e. A small insert in a larger repeat will have an overlap of 1. Range: 0-1" )
	parser.add_argument("--reportoverlaps", action="store_true", help="Report clusters of reference repeats which share hit locations overlapping >= 1bp with at least \ one other member of the cluster. Use as guide to curate and merge redundant reference repeats.")
	parser.add_argument("--writeoverlaps", action="store_true", help="If set with 'reportoverlaps', write groups of potentially redundant reference repeats \ to multi FASTA files for inspection.")
	#parser.add_argument("--doAlign",action="store_true",help="If set attempt to align reunited repeat families.")
	#parser.add_argument("--alnDir",type=str,default=None, help="Write aligned transposon families to this directory.")
	#parser.add_argument('--AlnTool',default="clustalW",choices=["clustalW","DIALIGN","MUSCLE"],help='Attempt to align clusters with this tool. \Note: Will first attempt alignment with DIALIGN if anchors set.')
	#parser.add_argument('--outAlnFormat',default="fasta",choices=["clustal","emboss","fasta","fasta-m10","ig","nexus","phylip","phylip-sequential","phylip-relaxed","stockholm"],help='Optional: Write alignment including reference sequence to file of format X.')
	return parser.parse_args()

def main():
	args = getArgs()
	# Check for required files + make output directories
	outDir = reunite.dochecks(args)
	# Import target genomes as {GenLabel:{ChromName:SeqRecord}}
	# Also return seq to genome map as {SeqName:GenLabel}
	genMaster,chr2Gen = reunite.importGenomes(args.genomes,args.genomeLabels)
	# Import reference repeat library, check for unique names.
	if args.blastTAB:
		refMaster = reunite.importRefseqs(args.repeats)
	elif args.rmOut:
		refMaster = reunite.importRefseqs(args.repeats,stripHash=True)
	# Import BLAST or RM hits, filter for query coverage and identity, determine orientation
	# Dict key by refName
	if args.blastTAB:
		allHits = reunite.importBLAST(infile=args.blastTAB, QuerySeq=refMaster,
								chr2Gen=chr2Gen, minCov=args.mincov,
								minID=args.minID, eVal=args.eVal)
	elif args.rmOut:
		allHits = reunite.importRM(infile=args.rmOut, QuerySeq=refMaster,
								chr2Gen=chr2Gen, minCov=args.mincov,
								minID=args.minID, score=args.minSW)
	else:
		sys.exit(1)
	# Find number of overlapping hits between any two refTEs
	if args.reportoverlaps:
		overlapDict = reunite.countIntersects(allHits,args.pOverlap)
		refClusters = reunite.groupOverlaps(overlapDict)
		reunite.writeOverlaps(outDir,refClusters)
		if args.writeoverlaps:
			reunite.writeRedundant(outDir,refClusters,refMaster)
	# Extract hits for each refRepeat, write clusters
	reunite.writeClusters(allHits,refMaster,genMaster,outDir,SkipZeros=args.onlyhits)
	# Write summary of hits per refRepeat found per target genome
	reunite.writeGenomeSummary(outDir, allHits, genMaster)
