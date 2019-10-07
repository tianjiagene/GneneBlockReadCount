#calculate gene  reads number, 
#for single end reads 
#usally using CDS region of refseq, also possible to calculate based on all exonic region
#use 04.cal_gene_Rnum.inc.pl to calculate
#can output a file for all junction reads mapping to gene information ($if_out_juncread_detail)

require "../SharedCodes/perl.fun.inc.pl";
use Getopt::Std; 

#######common setup
$debug=0; #1 0
$geno="mm9"; #hg18 mm8 rn4
$refsrc_geno=""; ###in most cases set to ""; when using refseq in another species set to the source species
$id_name="refseqid"; ###id name, if refseq, set "refseqid";  if reads cluster, rRNA, set "gene_id"; if ensembl, ensid
$if_directional_reads=1; ###!!!  1 (directional reads) or 0
$if_reverse_read_strand=0; #whether or not reverse read strand (eg. for reverse sequencing, read strand is always antisense to mRNA)
$if_reverse_read_strand_mate2=0; #whether or not reverse read strand for mate 2 of pair-end sequencing (eg. for reverse sequencing, read strand is always antisense to mRNA)
$countRnum_inRegion="CDS"; ### CDS (default); transcript ( include "CDS","3UTR","5UTR"); comprehensive
$ifout_detail_reads_info=0;

#######################samples:



getopt("sgilpeadrAuHDCEjBRo",\%args);
if ($args{s}){ #sample name string
	$study_name=$args{s}; ##mmMuscle2
	$geno=$args{g} if $args{g}; #
	$id_name=$args{i} if $args{i}; #
	$readslen=$args{l} ?$args{l}:0; #
	$countRnum_inRegion=$args{e} if $args{e}; #
	$sample_names_all=$args{a} if $args{a}; #
	$if_directional_reads=$args{d} if $args{d} ne ""; #
	$if_reverse_read_strand=$args{r} if $args{r} ne ""; #
	$if_reverse_read_strand_mate2=$args{R} if $args{R} ne ""; #
	@sample_names_arr=split(",| ",$sample_names_all) if $sample_names_all;
	$prefix=$args{p} ? $args{p} : ''; #process from which chromosome (in @all_chr_arr)
	$ano_file=$args{A} ? $args{A} : ''; 
	$ifout_detail_reads_info=$args{u} ? $args{u} : 0; 
	$r_GenoJuncMapNum_d=$args{D} ? $args{D} : "02.PE.ReadTable/"; #input read mapping info directory
	$contig_prefix=$args{C}?$args{C}:"";  
	$if_out_juncread_detail=$args{j} ne ""?$args{j}:"";
	$debug=$args{B} ne ""?$args{B}:0;
	$outdir=$args{o} ne ""?$args{o}:"$r_GenoJuncMapNum_d/../04.GeneReadNum/";
}

@r_GenoJuncMapNum_fexts=("GenoMapNum","JunReadNum");
if($args{E}){@r_GenoJuncMapNum_fexts=split(/\s/, $args{E}); }

################other setup
if(!$ano_file){
	if($id_name=~/refseq/ && !$refsrc_geno){
		$ano_file="../ReferenceDB/gene/02transcript_gene_ano/$geno.refflat.desc.txt"; #
	}elsif($id_name=~/ensid/){
		$ano_file="../ReferenceDB/gene/02transcript_gene_ano/$geno.ensGene.desc.txt";
	}
}
$cal_Rnum_cmd= "perl 04.cal_gene_Rnum.inc.pl -a \'$ano_file\'  -i $id_name -t $if_directional_reads -e $countRnum_inRegion -v $if_reverse_read_strand -R $if_reverse_read_strand_mate2 -j \"$if_out_juncread_detail\" ";
if($contig_prefix){
	$cal_Rnum_cmd="$cal_Rnum_cmd -C \"$contig_prefix\" ";
}

###################RUN
create_dir_ifNotExist("$outdir/");
foreach $samplename(@sample_names_arr){
	$r_GenoJuncMapNum_f="\'cat ". join(" ",map("$r_GenoJuncMapNum_d$samplename.$_",@r_GenoJuncMapNum_fexts)). " |\'";
	$ano_small_outf="$outdir/${prefix}ids.tbl";
	$readsnum_out_f="$outdir/${prefix}$samplename.ReadNum.tbl.temp";
	$detail_reads_info_outf=$ifout_detail_reads_info?"$outdir/${prefix}$samplename.detail.reads.txt":0 ;
	$cmd="$cal_Rnum_cmd -r $r_GenoJuncMapNum_f -o $readsnum_out_f -s $samplename -n $ano_small_outf -d $debug -f SE -l $readslen -u $detail_reads_info_outf ";
	print "cmd=$cmd\n";
	system ($cmd);
}

my $readsnum_allsample_out_f="$outdir/${prefix}ReadNum.tbl";
opendir(OUT_DIR, "$outdir/");
my @tmp_Rnum_files = sort grep(/$prefix.*\.ReadNum\.tbl\.temp$/,readdir(OUT_DIR));
unshift @tmp_Rnum_files, "${prefix}ids.tbl";
my $infiles_str=join(" ", map("$outdir/$_", @tmp_Rnum_files));
system("paste $infiles_str >$readsnum_allsample_out_f");
#delete temporary files
foreach my $tmp_file(@tmp_Rnum_files){
	#system("rm $outdir/$tmp_file");
}
