
##this pipeline (bash command lines) start with .bam files (i.e., output of RNA-seq read mapping using STAR program), map and count read to different regions of genes (exons, introns) 
#programs need to be installed in a linux server:
# perl, R
# STAR
# samtools

package_download_directory=/directory/to/your/downloaded/package/ #change this to the root directory of the package

study_name=example_project #project name
pipeline_root=$package_download_directory/Pipeline/ #this is root directory of the pipeline/scripts 
project_root=$package_download_directory/example_project #project folder (all the input and output files related to this project)
mapping_data_root=$project_root #use this directory to save large files (eg. fastq files and mapped bam files). expect a STAR_map folder with all .bam files under this directory


if_PE=0; #whether this is pair-end sequencing (0, no; 1, yes)
if_stranded=0; #whether the RNA library is stranded (contain gene strand information) (0, no; 1, yes)
if_read1_antiSense=0; #whether the first read is antisense to gene. This is usually true for TruSeq Stranded mRNA Library (second read is sense strand). (0, no; 1, yes)

##download hg19 genome fasta files from UCSC website:
##use bash commands in this file to download: $pipeline_root/ReferenceDB/ucsc/download_ucsc.sh
##save downloaded file in this folder: $pipeline_root/ReferenceDB/ucsc/genomes/hg19


####R.1, read mapping
##using STAR to map

# expect to have .bam files in folder $mapping_data_root/STAR_map
# if use STAR to map the reads, suggest use these settings: "STAR --runMode alignReads --outSAMattributes NH HI NM MD jM jI XS"
#note, in the example_project/STAR_map/ folder of this package, the bam files only have the reads mapped to chr2 (to reduce the file size)



####R.2, analyze RNA-seq data
####R2.1  find uniquely mapped reads with good quality from sam or bam file and output table format files for downstream analysis
study_name=example_project
samples="NC.41 CRC.41 NC.42 CRC.42"
cd $pipeline_root/RNAseq

p=0
for sample in $samples; do
  perl 01.Cal_geno_junc_readsnum_frSAM.pl -n $study_name -s "$sample"  -d "$mapping_data_root/STAR_map/" -o "$project_root/01.ReadTable/" \
    -e ".bam" -t STAR -q 10 -m 5 -u 0  &
  p=$(($p+1)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
done
wait

## "-q 10 -m 5" defines what kind of reads were used for downstream analysis
# -q  minimal MAPQ score required (default 10)
# -m  maximal number of mismatch per 100 nt of mapped length (default 5)



####R2.2 treat all reads as single end (SE) read, calculate junction read map to gene info (as well as gene expression based on CDS region of refseq transcripts)
##the output of this step will be used for splcing analysis for both PE and SE reads
geno=hg19
study_name=example_project
cd $pipeline_root/RNAseq
samples="NC.41 CRC.41 NC.42 CRC.42"
p=0
for sample in $samples; do
  perl 04.cal_gene_cds_readnum.pl -s $study_name -g $geno -i refseqid -p refseqcds_SE. -e CDS -d "$if_stranded" -r "$if_read1_antiSense" -a "$sample" \
    -D $project_root/01.ReadTable/ -j 1 -E "JunReadNum GenoMapNum" -o $project_root/04.GeneReadNum/ &
    p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait


#####Splicing analysis pipeline
####S.1 combine junction read mapping info 
geno=hg19
study_name=example_project
cd $pipeline_root/RNAseq_Splicing
samples="NC.41 CRC.41 NC.42 CRC.42"
Rscript 01.comb_junc_map_info.R -study_name $study_name -geno $geno -samples "$samples" -indir $project_root/04.GeneReadNum/refseqcds_SE. \
  -out_cbf $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
  -id_header refseqid
##this will output a file $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl which contains junction read count information for each sample


####S.2, cut genes to blocks (exon and intron part) based on Refseq annotation and novel splice sites identified from RNA-seq
geno=hg19
cd $pipeline_root/RNAseq_Splicing
perl 02.cut_gene3blocks_basedOnSplicing.pl  -g $geno -R 4 -S 2 -F 0.05 -i refseqid \
  -j $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
  -o $project_root/RNAseq_Splicing/02.geneBlocks

## this will create a file $project_root/RNAseq_Splicing/02.geneBlocks/blockBoundary.ano.tbl which contains gene block (exon, intron boundaries) information
# parameters "-R 4 -S 2 -F 0.05" are used to filter novel splice site (SS) identified from RNA-seq data using junction reads. 
#   -R  the minimal total read count (default 4), and
#   -F  the minimal relative fraction of read counts among the average of read counts in a gene (0.05 default) in at least two samples, and
#   -S  the minimal number of samples (default 2) have >0 read counts



####S.3, calculate gene block expression (exonic and intronic parts) 
geno=hg19
study_name=example_project
cd $pipeline_root/RNAseq
samples="NC.41 CRC.41 NC.42 CRC.42"
p=0
for sample in $samples; do
  perl 04.cal_gene_cds_readnum.pl -s $study_name -g $geno -i Gblock_id -p Gblock. -e transcript -d "$if_stranded" -r "$if_read1_antiSense" -a "$sample" \
    -D $project_root/01.ReadTable/ -j 0 -E "JunReadNum GenoMapNum" -o $project_root/RNAseq_Splicing/geneBlockReadNum/ \
    -A $project_root/RNAseq_Splicing/02.geneBlocks/geneBlock.flat.tbl  &
    p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait
## this will create a file $project_root/RNAseq_Splicing/geneBlockReadNum/Gblock.ReadNum.tbl which contains information of gene block supporting read counts for all samples

