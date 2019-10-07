
GeneBlockReadCount is a package of scripts for splicing analysis using RNA-seq data. The main function of this package is to count the read number for different regions (blocks) of a gene (eg. exons and introns). It uses not only Refseq-defined splice sites, but also novel splice sites from the RNA-seq data to define gene block boundaries. The output of this package can be used for further splicing analysis using DEXSeq or custome code etc.

GeneBlockReadCount use Perl and R scripts in a linux environment. 
Additional programs needed: samtools, STAR (used to create bam files as input of this package).
The current package only has gene structure (Refseq) annotation files for human genome under version hg19.

Find file example_project/example_GSE95132_CRC.sh for an example to use this pipeline. It contains the bash command lines to call different scripts in the pipeline. 
