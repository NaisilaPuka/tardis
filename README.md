tardis
======

Toolkit for Automated and Rapid DIscovery of Structural variants

Soylev, A., Kockan, C., Hormozdiari, F., & Alkan, C. (2017). Toolkit for automated and rapid discovery of structural variants. Methods, 129, 3-7. https://doi.org/10.1016/j.ymeth.2017.05.030

Soylev, A., Le, T., Amini, H., Alkan, C., & Hormozdiari, F. (2018). Discovery of tandem and interspersed segmental duplications using high throughput sequencing. bioRxiv, posted August 16, 2018. https://doi.org/10.1101/393694

Requirements
============

 * zlib   (http://www.zlib.net)
 * mrfast (https://github.com/BilkentCompGen/mrfast)
 * htslib (included as submodule; http://htslib.org/)
 * sonic  (included as submodule; https://github.com/calkan/sonic)

htslib also requires:

 * libbz2
 * liblzma


Fetching TARDIS
===============

	git clone https://github.com/BilkentCompGen/tardis.git --recursive

Compilation
===========

Type:

	make libs
	make
	cp tardis /path/to/your/favorite/binaries


SONIC file (annotations container)
==================================

SONIC files for some human genome reference versions are available at external repo: https://github.com/BilkentCompGen/sonic-prebuilt

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)
 	* Also download the reference genome at: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz. 
 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * ucsc_hg38.sonic: SONIC file for the human reference genome build 38.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.

Make sure that the same reference was used to align the reads beforehand (BAM file) and to create the SONIC file. The SONIC files and the reference FASTA files linked above are compatible.

Building the SONIC file
=======================

Please refer to the SONIC development repository: https://github.com/calkan/sonic/

However, you can still build the SONIC file using TARDIS:

	tardis --ref human_g1k_v37.fasta --make-sonic human_g1k_v37.sonic \
		--dups human_g1k_v37.segmental_duplications.bed \
		--gaps human_g1k_v37.assembly_gaps.bed \
		--reps human_g1k_v37.repeatmasker.out 

	

Running TARDIS - QUICK mode
===========================

	tardis -i myinput.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic  \
		--out myoutput


Running TARDIS - SENSITIVE mode (mrFAST Mappings)
=================================================

Sensitive mode uses mrFAST mappings (all possible mappings) with read-pair and read-depth signatures. 

	tardis -i myinput.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic  \
		--sensitive --out myoutput

This command first runs mrFAST and creates the DIVET file that contains all possible mappings of the reads in your BAM file. However, if you already have the DIVET files (there should be as many DIVET files as there are libraries in your BAM files), you can use --skip-mrfast

DIVET files should be inside the TARDIS directory or under the divet/ folder


All parameters
==============

	--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.
	--input [BAM files]        : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.
	--out   [output prefix]    : Prefix for the output file names.
	--ref   [reference genome] : Reference genome in FASTA format.
	--sonic [sonic file]       : SONIC file that contains assembly annotations.
	--read-count [int]         : # of clusters that a specific read can be involved in (Default is 10).
	--mei   ["Alu:L1:SVA"]     : List of mobile element names.
	--no-soft-clip             : Skip soft clip remapping.
	--no-interdup              : Skip interspersed duplication clustering.
	--resolved                 : Output sequence resolved vcf calls.
	--xa                       : Look for the alternative mapping locations in BWA.
	--first-chr [chr_index]	   : Start running from a specific chromosome [index in ref]
	--last-chr  [chr_index]	   : Run up to a specific chromosome [index in ref]

	Additional parameters for sensitive mode:

	--sensitive                : Sensitive mode that uses all map locations. Requires mrFAST remapping.
	--skip-mrfast              : Skip mrFAST mapping. Use this only if you already have the correct divet file. Sensitive mode only
	--threads                  : Number of threads for mrFAST to remap discordant reads.

	Additional parameters to build SONIC file within TARDIS:

	--make-sonic [sonic file]  : SONIC file that will contain the assembly annotations.
	--sonic-info [\"string\"]  : SONIC information string to be used as the reference genome name.
	--gaps  [gaps file]        : Assembly gap coordinates in BED3 format.
	--dups  [dups file]        : Segmental duplication coordinates in BED3 format.
	--reps  [reps file]        : RepeatMasker annotation coordinates in RepeatMasker format. See manual for details.
	
	Additional parameters for 10X Genomics Linked Reads (under development):

	--10x                      : Enable 10X Genomics Linked Reads mode.
	--output-hs                : Output the selected clusters homogeneity scores to the VCF file.

	Information:
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.


Converting output VCF file to BED
==============

	awk '! /\#/' out.vcf |\
	awk '{print $1"\t"($2-1)"\t"(substr($8,match($8,/SVLEN=[0-9]+/)+length("SVLEN="),RLENGTH-length("SVLEN="))+$2-1)}' > out.bed

Alternatively, use VCFlib: https://github.com/vcflib/vcflib

