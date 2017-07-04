# TE-Reunite

For a set of reference transposons, collect instances from one or more genomes and write as multi-FASTA files.
Repeat collections can be aligned and used to [build models capturing TE family diversity](http://hmmer.org/), or produce [deRIP'd reference sequences](https://github.com/Adamtaranto/deRIP2).

**Example usage:**

```bash

# Query reference repeats against a target genome
blastn -query RepeatLibrary.fa -evalue 0.01 -outfmt 6 -subject targetGenome.fa -out BLAST_Results.tab

# Launch Reunite on blast output file
# Report clusters of hits per reference repeat, and summary file of reference repeats with overlapping hit locations

./TE-Reunite.py -r RepeatLibrary.fa -b BLAST_Results.tab -c 0.8 -e 0.001 -i 0.7 -o output -g targetGenome.fa -l GEN_1 --reportoverlaps --writeoverlaps

```

## Options
  **-h, --help**  

  Show this help message and exit  

  **-v, --version**  

  Show program's version number and exit  

  **-r REPEATS, --repeats REPEATS**  

  Fasta formated library of reference repeats. Note: Names must be unique.  

  **-g GENOMES [GENOMES ...], --genomes GENOMES [GENOMES ...]**  

  Space delimited list of reference genomes used as
  BLAST targets. Fasta format. Note: Names of sequences
  must be unique across ALL genomes provided.  

  **-l GENOMELABELS [GENOMELABELS ...], --genomeLabels GENOMELABELS [GENOMELABELS ...]**  

  Space delimited list of genome labels. Note: Lables
  must be unique and same length as genome files list.  
  
  **-b BLASTTAB, --blastTAB BLASTTAB**  

  BLAST output format 6.  
  
  **--rmOut RMOUT**  

  Repeatmasker outfile.  
  
  **-s MINSW, --minSW MINSW**  
  
  Minimum Smith-Waterman alignment score if providing
  hits in RepeatMasker.out format.  
  
  **-c MINCOV, --mincov MINCOV**  

  Minimum coverage of reference sequence required for
  vaild hit. Range: 0-1  
  
  **-e EVAL, --eVal EVAL**  

  Minimum e-Value required for vaild hit.  
  
  **-i MINID, --minID MINID**
  Minimum sequence identity for vaild hit. Range: 0-1  

  **-o OUTDIR, --outDir OUTDIR**  

  Write reunited transposon families to this directory.  

  **--reportoverlaps**  

  Report clusters of reference repeats which share hit
  locations overlapping >= 1bp with at least one other
  member of the cluster. Use as guide to curate and
  merge redundant reference repeats.  

  **--writeoverlaps**  

  If set with 'reportoverlaps', write groups of potentially 
  redundant reference repeats to multi FASTA files for inspection.