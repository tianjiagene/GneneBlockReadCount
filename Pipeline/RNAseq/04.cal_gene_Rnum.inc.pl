

require "../SharedCodes/perl.fun.inc.pl";

use Getopt::Std; 

getopt("arflousnditevHCjR",\%args);
$geneAno_f=$args{a}; #
$reads_mapInfo_f=$args{r}; # could use "cat file1 file2 |" to combine two or more files
$reads_mapInfo_f_format=$args{f}; #"PE.cov" "SE" "SE"
$readslen=$args{l}; #reads length (only for SE.geno)
$readsnum_out_f=$args{o}; # gene reads num based on CDS region
$detail_reads_info_outf=$args{u}; # other reads mapping info: intron, UTR5e-, intergenic reads
$samplename=$args{s}; #
$geneAno_outf=$args{n}; #short gene anotation out file (used for combining to final out file)
$debug=$args{d}; #if debug
$id_name=$args{i}; #id name, if refseq, set "refseqid"; if reads cluster, set "gene_id"
$if_directional_reads=$args{t}; 
$countRnum_inRegion=$args{e}; #CDS (default); transcript ( include "CDS","3UTR","5UTR"); comprehensive (antisense_exon, antisense_intron, sense_exon, sense_intron)
$if_reverse_read_strand=$args{v}; #whether or not reverse read strand for single end sequencing or mate 1 of pair end sequencing (eg. for reverse sequencing, read strand is always antisense to mRNA)
$if_reverse_read_strand_mate2=$args{R}; #whether or not reverse read strand for mate 2 of pair-end sequencing (eg. for reverse sequencing, read strand is always antisense to mRNA)
$if_out_juncread_detail=$args{j} ne ""?$args{j}:"";



$id_name="refseqid" if !$id_name; #default: treat refseq file
$countRnum_inRegion="CDS" if !$countRnum_inRegion; # CDS (default); transcript ( include "CDS","3UTR","5UTR")
$contig_prefix=$args{C}?$args{C}:"";  

print "####\n";
foreach $var(split(/,/,"geneAno_f,reads_mapInfo_f,reads_mapInfo_f_format,readslen,readsnum_out_f,detail_reads_info_outf,samplename,geneAno_outf,if_directional_reads,if_reverse_read_strand,debug,countRnum_inRegion,if_out_juncread_detail")){
	print "$var=$$var\n";
}

#####setup
$reads_mapInfo_f=~s/\'//g;

$logfile="$readsnum_out_f.log";
$chr_block_size=5000;
$ifuse_1endMapped_reads=1; #whether or not use only one end mapped reads
$utr_ext_len=4000; 

%readsTypes2sco_h=( #higher score represent higher priority (applicable to RNA-seq); + - represent strand relative to reference gene only useful for directional reads 
	"CDS_S"=>18,
	"CDS"=>17,
	"3UTR_S"=>16,
	"3UTR"=>15,
	"5UTR_S"=>14,
	"5UTR"=>13,
	"exon_S"=>12.9,
	"exon"=>12.8,
	"intron_S"=>12,
	"intron"=>11,	
	"UTR3e_S"=>10,
	"UTR3e"=>9,
	"intron_A"=>8,
	"UTR5e_A"=>7,
	"UTR5e"=>6,
	"UTR5e_S"=>5,
	"UTR3e_A"=>4,
	"3UTR_A"=>3,
	"CDS_A"=>2,
	"exon_A"=>1.9,
	"5UTR_A"=>1,
	"intergenic"=>0,
	""=> -1
);
$detail_reads_outtypes=".*";
@comprehensive_out_regions=("CDS_S","5UTR_S","3UTR_S","intron_S","UTR3e_S","CDS_A","5UTR_A","3UTR_A","intron_A","UTR3e_A","UTR5e_A");


###################RUN
create_dir_ifNotExist($readsnum_out_f);
open (LOGF,">$logfile") || die "error $logfile\n";



#1, load gene anotation file into hash %chro2reginfo_hash
open (GENE_ANO, $geneAno_f) || die "error $geneAno_f\n";
open (GENE_OUT, ">$geneAno_outf")  || die "error $geneAno_outf\n";
print " open $geneAno_f\n write $geneAno_outf\nheaders:\n";
print LOGF " open $geneAno_f\n write $geneAno_outf\n";
print GENE_OUT "$id_name	cds_len	transcript_len\n";
my $rowi=0;
while(<GENE_ANO>){
	$rowi++;
	chomp;
	if(/^gene_symbol|$id_name	/i){ #if this changed, also change another line below in step 3. # gene_symbol refseqid contig  strand  transc_start transc_end cds_start cds_end  exon_num exon_starts exon_ends gene_id alias   gene_desc
		@headername_arr=split(/\t/);
		print $_."\n";
	}else{
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}


		#modify last exon end position
		if($gene3end_redefine_f){
			($transc_start,$transc_end)=modify_exon3end($contig, $strand, $transc_start,$transc_end, $exon_starts,$exon_ends);
		}
		
		if($cds_start==$cds_end || !$cds_start){
			#($cds_start,$cds_end)=($transc_start,$transc_end);
			$cds_start=0; $cds_end=0; #changed 3/17/2014
		}
		($cds_len,$transcript_len)=cal_refseq_genelen($cds_start,$cds_end,$exon_starts,$exon_ends);
		$utr5e_pos=$transc_start-$utr_ext_len;
		$utr3e_pos=$transc_end+$utr_ext_len;
		if($strand eq "-" || $strand eq "-1"){
			($cds_start,$cds_end, $transc_start,$transc_end, $utr5e_pos,$utr3e_pos)=($cds_end,$cds_start, $transc_end,$transc_start, $utr3e_pos,$utr5e_pos ); ####for refseq###
		}
		
		print GENE_OUT "$$id_name	$cds_len	$transcript_len\n";
		if($contig_prefix){$contig="$contig_prefix$contig";} #used for modified contig name eg. yeast genome in mouse and yeast combined genome index
		

		if($cds_start){
			create_chr_reg_info_hash($cds_start,$cds_end, "CDS.$contig", "$exon_starts|$exon_ends|$strand", $chr_block_size, join("	",(sort{$a<=>$b}($cds_start,$cds_end)))."	$$id_name	$rowi" );
			create_chr_reg_info_hash($transc_start,$cds_start, "5UTR.$contig", "$exon_starts|$exon_ends|$strand", $chr_block_size, join("	",(sort{$a<=>$b}($transc_start,$cds_start)))."	$$id_name	$rowi" );
			create_chr_reg_info_hash($cds_end,$transc_end, "3UTR.$contig", "$exon_starts|$exon_ends|$strand", $chr_block_size, join("	",(sort{$a<=>$b}($cds_end,$transc_end)))."	$$id_name	$rowi" );
			create_chr_reg_info_hash($utr5e_pos, $transc_start, "UTR5e.$contig",  "$utr5e_pos|$transc_start|$strand",  $chr_block_size, join("	",(sort{$a<=>$b}($utr5e_pos, $transc_start)))."	$$id_name	$rowi" );
			create_chr_reg_info_hash($transc_end, $utr3e_pos, "UTR3e.$contig",  "$transc_end|$utr3e_pos|$strand",  $chr_block_size, join("	",(sort{$a<=>$b}($transc_end, $utr3e_pos)))."	$$id_name	$rowi" );
		}else{
			create_chr_reg_info_hash($transc_start,$transc_end, "exon.$contig", "$exon_starts|$exon_ends|$strand", $chr_block_size, join("	",(sort{$a<=>$b}($transc_start,$transc_end)))."	$$id_name	$rowi" );
		}
		$cds_start="";
		$cds_end="";
	}	
}
close GENE_ANO;
close GENE_OUT;


#2, read reads map info and count number mapped to gene cds region to %gene2readsnum_hash
my %gene2readsnum_hash=();
print " open $reads_mapInfo_f\n";
open (READS_COV, $reads_mapInfo_f) || die "error reads_mapInfo_f=$reads_mapInfo_f\n";
print LOGF " open $reads_mapInfo_f\n";
if($detail_reads_info_outf){
	open (OTH_READSOUT, ">$detail_reads_info_outf") || die "error detail_reads_info_outf=$detail_reads_info_outf\n";
	print OTH_READSOUT "chromosome	reads_strand	posfr	posto	readsnum	readsType	geneInfo\n";
	print LOGF " write $detail_reads_info_outf\n";
}
my $rowj=0;
my $ifJuncReadFile=0; my %juncReadInfo_h=(); my $junc_strand="";
while(<READS_COV>){
	chomp;
	if(/^chromosome|readsnum/i){ 
		#PE.cov file:  chromosome  posfr   posto   readsnum
		#reads.cov.tbl chromosome      mate1strand     posfr   posto   readsnum
		#SE.geno file: rname   reads_strand    posi    readsnum
		# rname   reads_strand    posfr   posto   readsnum
		#SE.junc file: rname   junc_strand  juncpos5  juncpos3  readsnum  exonreg5lens
		s/rname/chromosome/;
		if(/juncpos[53]/){$ifJuncReadFile=1;}else{$ifJuncReadFile=0;}
		s/posi|juncpos5/posfr/;
		s/juncpos3/posto/;
		s/junc_strand|mate1strand/reads_strand/;
		@headername_arr2=split(/\t/);
	}else{
		$rowj++;
		@temp_arr = split(/\t/); #chr1    3522250 3522311 2
		$posto=0;
		for($i=0;$i<@headername_arr2;$i++){
			${$headername_arr2[$i]}=$temp_arr[$i];
		}
		if($ifJuncReadFile){$junc_strand=$reads_strand;}
		if(!$if_directional_reads){$reads_strand="";}
		if($if_directional_reads && $if_reverse_read_strand && !$ifJuncReadFile){
			$reads_strand=($reads_strand=~/\-/ ? '+' : '-');
		}
		next if ($posfr==0 || $chromosome eq "*"); #unmapped in eland output file
		

		$posto=$posfr+$readslen-1 if (!$posto && $readslen>0);
		die "read map end not defined in file or read length not defined!\n" if !$posto;
		if (!$ifuse_1endMapped_reads && $readsnum eq "1" && $reads_mapInfo_f_format eq "PE.cov"){
			$counts{"reads.unused_1endmapped"}+=$readsnum;
			next;
		}
		if($rowj>100000 && $debug){
			print "debug skipped at row $rowj of $reads_mapInfo_f\n";
			last;
		}; 
		if($rowj % 1000000==0 || eof){
			print  return_time()."|read reads_mapInfo_f file $samplename line $rowj\n";
			print  LOGF return_time()."|read reads_mapInfo_f file $samplename line $rowj\n";
		}
		


		my $posfrom_block_num=int($posfr/$chr_block_size);
		my $posto_block_num=int($posto/$chr_block_size);
		($posfrom_block_num,$posto_block_num)=sort {$a<=>$b} ($posfrom_block_num,$posto_block_num);


		my %treated_ids=(); # one reads can only contribute one refseq id once
		my $readsType_final="intergenic"; 
		my $transcript2readsType_now=();
		my $geneInfo_now=""; 
		my $geneInfo_final="";
		
		foreach $reg_type(("CDS","3UTR","5UTR","exon","UTR3e","UTR5e")){
			foreach my $block_num(($posfrom_block_num..$posto_block_num)){
				foreach my $key2 ( keys %{$chro2reginfo_hash{"$reg_type.$chromosome:$block_num"}} ){
					($reg_start,$reg_end,$refseqid,$rowi)=split(/	/,$key2);
					
					next if ($treated_ids{$rowi}); #one read only contribute one gene once
					$geneinfo_str=$chro2reginfo_hash{"$reg_type.$chromosome:$block_num"}{$key2};
					($exon_starts,$exon_ends,$gene_strand)=split(/\|/,$geneinfo_str);
					@exon_start_arr=split(/,/,$exon_starts);
					@exon_end_arr=split(/,/,$exon_ends);

					next if ($ifJuncReadFile  && $junc_strand ne $gene_strand); #junction read always have direction
					$reads2cds = judge_2region_relation($reg_start,$reg_end, $posfr,$posto, $gene_strand);
					next if ($reads2cds!~/0|ov-1|=/); #rough judgement
					$rel_strand=($reads_strand? ( ($gene_strand eq $reads_strand)?"_S":"_A" ) : ""); 
					
					#map to at least intron
					$transcript2readsType_now{$refseqid}="intron$rel_strand"; #read in feature bundary, but not sure map to exon, assign to intron first
					$geneInfo_now=$refseqid;

					for(my $i=0;$i<@exon_start_arr;$i++){ #check every exons, rule out reads in the intron region only
						$exon_start=$exon_start_arr[$i];
						$exon_end=$exon_end_arr[$i];
						$reads2exon = judge_2region_relation($exon_start,$exon_end, $posfr,$posto, $gene_strand);
						if($reads2exon=~/0|ov-1|=/){ #reads inside or overlap with exon (not allow read to overlap downstream of feature)
							$treated_ids{$rowi}=1; #if overlap with one gene, no need to search this gene any more
							$transcript2readsType_now{$refseqid}="$reg_type$rel_strand";
							$geneInfo_now=$refseqid;
							last; #no need to search other exons for the same gene
						}
					}
					
					if($readsTypes2sco_h{$transcript2readsType_now{$refseqid}}>$readsTypes2sco_h{$readsType_final}){
						$readsType_final=$transcript2readsType_now{$refseqid}; $geneInfo_final=$geneInfo_now;
					}
					
					#record read map to transcript and read type
					if($transcript2readsType_now{$refseqid} ne "intergenic"){
						$gene2readsnum_hash{$transcript2readsType_now{$refseqid}."	$refseqid	$rowi"}+=($ifuse_1endMapped_reads?$readsnum:1); #one read may count several times for several genes, but only count once for one gene
					}

					##record junction read information
					if($if_out_juncread_detail && $ifJuncReadFile  && $transcript2readsType_now{$refseqid}=~/CDS|[35]UTR|exon|intron|UTR3e|UTR5e/ && $junc_strand eq $gene_strand){
						$juncReadInfo_h{$transcript2readsType_now{$refseqid}."	$refseqid	$rowi"}{"$posfr	$posto"}+=$readsnum; #ss5 and ss3 (exonic 1-based coordinates)
					}
	
				}#end foreach $key2
			}#end foreach	$block_num	 	
		}#end foreach $reg_type
		
		my $count_key="reads.$readsType_final";
		if($if_out_juncread_detail){
			$count_key="reads.".($ifJuncReadFile?"Junc.":"nonJunc.")."$readsType_final";
		}
		$counts{$count_key}+=$readsnum; 
		if($detail_reads_info_outf && $readsType_final=~/$detail_reads_outtypes/){ #output OTH_READSOUT
			print OTH_READSOUT "$chromosome	$reads_strand	$posfr	$posto	$readsnum	$readsType_final	$geneInfo_final\n";
		}
	}#end else (finish sam file one line)	
}
close READS_COV;
close OTH_READSOUT if($detail_reads_info_outf);


#3, output reads number
open (GENE_ANO, $geneAno_f) || die "error $geneAno_f\n";
open (OUT, ">$readsnum_out_f")  || die "error $readsnum_out_f\n";
print " open $geneAno_f\n write $readsnum_out_f\n";
print LOGF " open $geneAno_f\n write $readsnum_out_f\n";
if($countRnum_inRegion eq "comprehensive"){	
	print OUT join("	", map("num_${_}_$samplename",@comprehensive_out_regions))."\n";
}else{
	print OUT "num_$samplename\n";
}

$rowi=0;
while(<GENE_ANO>){
	$rowi++;
	chomp;
	if(/^gene_symbol|$id_name	/i){ #gene_symbol refseqid contig  strand  transc_start transc_end cds_start cds_end  exon_num exon_starts exon_ends gene_id alias   gene_desc
		@headername_arr=split(/\t/);
	}else{
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}
		my $gene_Rnum=0;
		if($countRnum_inRegion eq "CDS"){ #default
			$gene_Rnum=$gene2readsnum_hash{"CDS	$$id_name	$rowi"}+$gene2readsnum_hash{"CDS_S	$$id_name	$rowi"} + 
			  $gene2readsnum_hash{"exon	$$id_name	$rowi"}+$gene2readsnum_hash{"exon_S	$$id_name	$rowi"};
		}elsif($countRnum_inRegion eq "comprehensive"){
			$gene_Rnum=join("	", map($gene2readsnum_hash{"$_	$$id_name	$rowi"}, @comprehensive_out_regions));
		}else{ #transcript
			$gene_Rnum=$gene2readsnum_hash{"CDS	$$id_name	$rowi"}+$gene2readsnum_hash{"3UTR	$$id_name	$rowi"}+$gene2readsnum_hash{"5UTR	$$id_name	$rowi"}+
				$gene2readsnum_hash{"CDS_S	$$id_name	$rowi"}+$gene2readsnum_hash{"3UTR_S	$$id_name	$rowi"}+$gene2readsnum_hash{"5UTR_S	$$id_name	$rowi"}+
				$gene2readsnum_hash{"exon	$$id_name	$rowi"}+$gene2readsnum_hash{"exon_S	$$id_name	$rowi"};
		}
		print OUT "$gene_Rnum\n";
	}	
}
close GENE_ANO;
close OUT;

###output junction read info
if($if_out_juncread_detail){
	my $out_junc_info_f="$readsnum_out_f.junc2gene.tbl";
	open (JUNC_OUT, ">$out_junc_info_f") || die "error write $out_junc_info_f\n";
	print JUNC_OUT "readType	$id_name	row	juncpos5	juncpos3	num_$samplename\n";
	foreach my $gene_info_str(sort keys %juncReadInfo_h){
		foreach my $junc_coor(sort keys %{$juncReadInfo_h{$gene_info_str}} ){
			print JUNC_OUT "$gene_info_str	$junc_coor	".$juncReadInfo_h{$gene_info_str}{$junc_coor}."\n";
		}

	}
	close JUNC_OUT;
}

#numbers
print "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
print LOGF "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
close LOGF;


sub judge_moreProper_Readstype{
	($type1,$type2)=@_;
	$readsTypes2sco_h{$type1} >= $readsTypes2sco_h{$type2}?$type1:$type2 ;
}

sub modify_exon3end{ # input are all ucsc raw data, $transc_start always<$transc_end, start position is 0 based
	my($contig, $strand, $transc_start,$transc_end, $exon_starts,$exon_ends)=@_;
	my @exon_start_arr=split(/,/,$exon_starts);
	my @exon_end_arr=split(/,/,$exon_ends);	
	
	if($strand eq "-" || $strand eq "-1"){
		if(my $redefined_pos3=$redefine3end_h{$contig}{$strand}{ $exon_end_arr[1] }){
			$transc_start=$redefined_pos3;
			$counts{"num.gene3end_redefined.actural"}++;
		}
	}else{ # + strand
		if( my $redefined_pos3=$redefine3end_h{$contig}{$strand}{ ( $exon_start_arr[-1]+1) } ){
			$transc_end=$redefined_pos3;
			$counts{"num.gene3end_redefined.actural"}++;
		}
	}
	return ($transc_start,$transc_end);
}
