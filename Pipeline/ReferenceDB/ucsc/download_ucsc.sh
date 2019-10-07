
echo $pipeline_root


######### download the genome sequence files (fasta format) (these files only need to download once)
for geno in hg19; do   # hg19 mm9
 echo __________ download the genome sequence of $geno __________
 dir=$pipeline_root/ReferenceDB/ucsc/genomes/$geno
 mkdir -p $dir
 cd $dir
 rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/$geno/chromosomes/ ./
done

#unzip and combine fasta files:
for geno in hg19; do   # hg19 mm9
 echo __________ unzip and combine genome sequence for $geno __________
 dir=$pipeline_root/ReferenceDB/ucsc/genomes/$geno
 cd $dir
 gunzip *gz
 cat chr*.fa >$geno.fa
done

