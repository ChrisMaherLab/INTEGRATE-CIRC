# INTEGRATE-Circ
***
INTEGRATE-Circ is a fusion junction detection tool, capable of identifying gene fusion transcripts and various isoforms that may result from these fusions, including fusion-derived circRNA (fcircRNA). INTEGRATE-Circ is an extension of the algorithm employed by INTEGRATE, which can be found [here](https://github.com/ChrisMaherLab/INTEGRATE).

INTEGRATE-Circ is developed at the [Christopher Maher Lab](http://www.maherlab.com/) at [Washington University in St. Louis](http://www.wustl.edu).
***

## Installation
First, clone the INTEGRATE-Circ repository:
```
cd PATH_TO_TOOL
git clone https://github.com/ChrisMaherLab/INTEGRATE-CIRC.git
```
Next, compile the tool (requires cmake version >= 2.8):
```
cd INTEGRATE-Circ
mkdir build
cd build
cmake ..
make
```
An executable file should now be located at ``` PATH_TO_TOOL/INTEGRATE-CIRC/build/bin/Integrate-Circ ```. To make this file easily accessible from other locations, update the $PATH variable using ```export PATH="PATH_TO_TOOL/INTEGRATE-CIRC/build/bin:$PATH" ``` at the command line.

## Preparing Input
Besides user-provided sequencing data, two input files are required for INTEGRATE-Circ.

1. Reference genome fasta file
2. Ensembl annotation file

Examples of each can be found in the example-data directory.

## Quickstart
Prior to analyzing sequencing data, INTEGRATE-Circ must perform a BWT (Burrows-Wheeler transform) of the reference genome. This operation only needs to be performed once and can be done using the following command. If no output directory is specificed, output defaults to ```./bwts```.
```
Integrate-Circ mkbwt /path/to/reference.fa -dir /path/to/bwt_directory
```

After running the above command, sequencing analysis can be performed:
```
Integrate-Circ fusion \
	/path/to/reference.fa \
	/path/to/annotation_file.txt \
	/path/to/bwt_directory \
	/path/to/mapped_RNA_reads.bam \
	/path/to/unmapped_RNA_reads.bam \
	/path/to/mapped_tumor_DNA.bam \
	/path/to/mapped_normal_DNA.bam
```
By default, all outputs are placed in the current working directory.

If mapped and unmapped RNA reads are in the same file, simply provide that same file twice. Both DNA files are optional. ```Integrate-Circ fusion -h``` can be run for a more complete description of parameters, which are also described below.

## Output
The INTEGRATE-circ fusion command outputs the following files:
1. summary.tsv: Contains a summary of all called fusions as well as any associated splice variants, including backsplices
2. exons.tsv: Contains information regarding fusion junctions, useful for creating primers for validation sequencing
3. breapoints.tsv: Describes fusion breakpoints as determined by RNA and DNA data
4. bk_sv.vcf: Describes fusion breakpoints in vcf format
5. fusions.bedpe: Describes fusion and splice junctions in SMC-RNA bedpe format
6. reads.txt: Describes information about all sequencing reads that support the identified junctions
7. fcirc.txt: Descibes fcircRNAs. Columns 1 and 2 provide ID and fusion gene information. Column 3 describes the backsplice acceptor, column 4 describes the backsplice donor, column 5 descibes the 5' fusion junction and column 6 describes the 3' fusion junction. Note that in a geneA::geneB fusion, columns 3 and 5 describe locations in geneA and columns 4 and 6 describes locations in geneB.

## Example
Analysis of the example-data can be performed as follows. The example data includes simulated reads consistent with a TMPRSS2::ERG gene fusion that also creates an fcircRNA.
```
cd PATH_TO_TOOL
mkdir example-data/bwts
Integrate-Circ mkbwt example-data/chr21.fa -dir example-data/bwts
Integrate-Circ fusion \
	example-data/hg19.chr21.fa \
	example-data/hg19.chr21.annotation.txt \
	example-data/bwts \
	example-data/RNA.Example.AllReads.bam \
	example-data/RNA.Example.AllReads.bam
```
By default, all outputs are placed in the current working directory. For visualization of the TMPRSS2::ERG fcircRNA, the output fcirc.txt file can be supplied as input to the most recent version of the [INTEGRATE-vis tool developed by the Chris Maher Lab](https://github.com/ChrisMaherLab/INTEGRATE-Vis).

## Additional parameters
Running ```INTEGRATE-Circ mkbwt``` (without parameters) gives the following parameter options:
```
    Integrate-Circ mkbwt (options) reference.fasta

    options:

            -mb  integer  :     sequences in the reference fasta that are shorter than this value        default: 10000000
                                are not included in the evaluation of repetitive reads.   
            -dir string   :     directory to store the BWTs.                                             default: ./bwts

```
Running ```INTEGRATE-Circ fusion -h``` gives the following information:
```
Integrate-Circ fusion (options) reference.fasta annotation.txt directory_to_bwt accepted_hits.bam unmapped.bam (dna.tumor.bam dna.normal.bam)

options: -cfn      integer : Cutoff of spanning RNA-Seq reads for fusions with non-canonical
                             exonic boundaries.                                                         default: 3
         -rt       float   : Normal dna / tumor dna ratio. If the ratio is less than
                             this value, then dna reads from the normal dna data set 
                             supporting a fusion candidates are ignored.                                default: 0.0
         -minIntra integer : If only having RNA reads, a chimera with two adjacent
                             genes in order is annotated as intra_chromosomal rather than 
                             read_through if the distance between the two genes is larger than
                             this value.                                                                default: 400000
         -minW     float   : Mininum weight for the encompassing rna reads on an edge.                  default: 2.0
         -mb       integer : See subcommand "mkbwt".
                             This value can be larger than used by mkbwt.                               default: 10000000
         -minDel   int     : minimum size of a deletion that can cause a fusion.                        default: 5000
         -reads    string  : File to store all the reads.                                               default: reads.txt
         -sum      string  : File to store summary.                                                     default: summary.tsv
         -ex       string  : File to store exons for fusions with canonical exonic boundaries.          default: exons.tsv
         -bk       string  : File to store breakpoints                                                  default: breakpoints.tsv
         -vcf      string  : File to store breakpoints in vcf format                                    default: bk_sv.vcf
         -fcirc    string  : File to store fcirc results in                                             default: fcirc.txt
         -bedpe    string  : File to store all fusions in SMC-RNA bedpe format                          default: junctions.bedpe
         -bacc     integer : max difference between spanning reads and annotation to decide canonical.  default: 1
         -largeNum integer : if a gene shows greater or equal to this number, remove it from results.   default: 4
         -sample   string  : sample name                                                                default: sample

This version of Integrate-Circ works in the following situations:
(1)having rna tumor, dna tumor, dna normal
(2)having rna tumor, dna tumor
(3)having rna tumor

Integrate-Circ will only use sequences in reference.fasta. 
Chr names with and without "chr" are regarded as the same, e.g. chr1 = 1.
The rna and dna bams can be from alignments mapped to different reference files with different order of the sequences and their names with or without "chr". However, The versions should be the same, e.g. hg19. (Also, the same as in annotation.)
The tumor and normal dna bams should be mapped to the same reference file.

For rna tumor: accepted_hits.bam is a bam file containing mapped rna reads. unmapped.bam is a bam contains the not mapped rna reads. If they have been merged into one bam, just use merged.bam twice in the command line.

For dna bams: If solt-clips are provided, then Integrate-Circ is trying to search rearrangement breakpoints, otherwise, only paired reads may be included in the analysis.

If having rna normal only or having both rna and dna normal data sets. These data sets can be run to find non somatic events.
e.g. Integrate-Circ fusion -normal (options) reference.fasta annotation.txt directory_to_bwt accepted_hits.normal.bam unmapped.normal.bam (dna.normal.bam)

## Expected run time 
In our initial benchmarking, INTEGRATE-Circ  was run on a big memory blade with 32 Intel Xeon CPU E5-2640s with 400G of memory and was able to process input whole genome sequencing data with ~200 million reads in approximately 1.5 hours. This is comparable to the original INTEGRATE tool, which was also run on 32 Intel Xeon CPU E5-2640s with 400G of memory and was able to analyze ~350 million sequencing reads in approximately 8 hours.

```
