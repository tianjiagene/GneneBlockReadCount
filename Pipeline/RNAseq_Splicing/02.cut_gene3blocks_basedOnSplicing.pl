##cut a gene to blocks based on splicing information from RNA-seq data.

require "../SharedCodes/perl.fun.inc.pl";
use Getopt::Std; 

$geno="hg19"; #hg19 mm9 rn4
$id_name="refseqid"; ###id name, if refseq, set "refseqid";  
@AS_types=("a5ss","a3ss","se");

%asType2_anoHeaders_h=(
	a5ss=>"exon_i	exon_uc	exon_dc	junc_i	junc_e	exon_i_len	exon_i_inFrame	exon_i_seq	exon_i_aa	exon_i_stopCodon	allJunc_ReadNum",
	a3ss=>"exon_i	exon_uc	exon_dc	junc_i	junc_e	exon_i_len	exon_i_inFrame	exon_i_seq	exon_i_aa	exon_i_stopCodon	allJunc_ReadNum",
	se=>"exon_i	exon_uc	exon_dc	junc_ui	junc_di	junc_e	exon_i_len	exon_i_inFrame	exon_i_seq	exon_i_aa	exon_i_stopCodon	allJunc_ReadNum"
);

getopt("sgiAjoSRFNPa",\%args);
$geno=$args{g} if $args{g}; #
$contig_file_dir="../ReferenceDB/ucsc/genomes/$geno";
$contig_f="";
if(-e "$contig_file_dir/$geno.fa"){ $contig_f="$contig_file_dir/$geno.fa"; }
$id_name=$args{i} ? $args{i} : "refseqid"; #
$geneAno_f=$args{A} ? $args{A} : '';
$junction_inf=$args{j} ? $args{j} : "01.comb_junc_map_info/combine.junc2gene.tbl";
$out_root=$args{o} ? $args{o} : "02.geneBlocks";
$outGeneBlock_flatf="$out_root/geneBlock.flat.tbl";
$outGeneBlock_anof="$outGeneBlock_flatf.ano";
$outAllJunc_f="$out_root/AllJunc.tbl";
$outBlockBoundary_f="$out_root/blockBoundary.ano.tbl";
$outAS_anofroot=$args{a} ne "" ? $args{a} : "$out_root/ASano.";
$novelJunc_SuppSampleMin=$args{S} ne "" ? $args{S} :2;
$novelJunc_SuppReadMin=$args{R} ne "" ? $args{R} :3;
$novelJunc_SuppRelReadMin=$args{F} ne "" ? $args{F} :0.05; #F means a Fraction
$sampleNames_str=$args{N} ne "" ? $args{N} : "";
$sampleNames_autoDetect_patt=$args{P} ne "" ? $args{P} : "^num_";

if(!$geneAno_f){
	if($id_name=~/refseq/ && !$refsrc_geno){
		$geneAno_f="../ReferenceDB/gene/02transcript_gene_ano/$geno.refflat.desc.txt"; #
	}elsif($id_name=~/ensid/){
		$geneAno_f="../ReferenceDB/gene/02transcript_gene_ano/$geno.ensGene.desc.txt";
	}else{
		die "error: geneAno_f not defined!\n";
	}
}


$log_f="$outGeneBlock_flatf.log";
create_dir_ifNotExist($log_f);



####RUN
open (LOGF, ">$log_f") || die "error write $log_f\n";

##1, read geneAno_f and load all known splicing and gene information: %geneCut_point_h; %kn_ss_h;  %all_juncs; %gene2junc_h
open (GENE_ANO, $geneAno_f) || die "error read geneAno_f=$geneAno_f\n";
print " open $geneAno_f\n headers:\n";
print LOGF " open $geneAno_f\n";
my %all_juncs=();
my $rowi=0;
while(<GENE_ANO>){
	$rowi++;
	chomp;
	if(/^gene_symbol|$id_name	/i){ #eg: gene_symbol refseqid contig  strand  transc_start transc_end cds_start cds_end  exon_num exon_starts exon_ends gene_id alias   gene_desc
		@headername_arr=split(/\t/);
		print $_."\n";
	}else{
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}
		
		if($cds_start==$cds_end || !$cds_start){
			$cds_start=0; $cds_end=0; #changed 3/17/2014
		}else{
			$cds_start++; #convert to 1-based coordinates
		}
		$transc_start++; ##convert to 1-based coordinates
		if($strand eq "-" || $strand eq "-1"){
			($cds_start,$cds_end, $transc_start,$transc_end, $utr5e_pos,$utr3e_pos)=($cds_end,$cds_start, $transc_end,$transc_start, $utr3e_pos,$utr5e_pos ); 
		}
		if($cds_start){
			$geneCut_point_h{"$contig	$strand	$$id_name	$rowi"}{"$cds_start"}="cs	cs";
			$geneCut_point_h{"$contig	$strand	$$id_name	$rowi"}{"$cds_end"}="ce	ce";
		}
		#find known ss, load %kn_ss5_h; %kn_ss3_h
		my @exon_start_arr=split(/,/,$exon_starts);
		my @exon_end_arr=split(/,/,$exon_ends);	
		map($exon_start_arr[$_]++,  0..(scalar @exon_start_arr-1) ); #convert to 1-based coordinates
		if($strand eq "-" || $strand eq "-1"){
			my @tmpArr=reverse(@exon_start_arr);
			@exon_start_arr=reverse(@exon_end_arr);
			@exon_end_arr=@tmpArr;
		}
		#print "$$id_name $contig:$strand; cds_len=$cds_len,transcript_len=$transcript_len\n\$exon_starts=$exon_starts\n\$exon_ends=$exon_ends\n";
		foreach my $exon_i(1..@exon_start_arr){
			my $ss3=$exon_start_arr[$exon_i-1];
			my $ss5=$exon_end_arr[$exon_i-1];
			if ($exon_i>1){
				my $ups_ss5=$exon_end_arr[$exon_i-2];
				$all_juncs{"$contig	$strand	$ups_ss5	$ss3"}++; #all 1-based exonic coordinates
				$gene2junc_h{"$contig	$strand	$$id_name	$rowi"}{"$ups_ss5	$ss3"}++;
			}
			
			$geneCut_point_h{"$contig	$strand	$$id_name	$rowi"}{"$ss3"}=$exon_i>1 ? "ss3	E$exon_i" : "ts	E$exon_i";
			$geneCut_point_h{"$contig	$strand	$$id_name	$rowi"}{"$ss5"}=$exon_i<@exon_start_arr ? "ss5	E$exon_i" : "te	E$exon_i";
			$kn_ss_h{'3'}{"$contig:$strand:$ss3"}{$$id_name}="E$exon_i" if $exon_i>1;
			$kn_ss_h{'5'}{"$contig:$strand:$ss5"}{$$id_name}="E$exon_i" if $exon_i<@exon_start_arr;
		}
		if($cds_start && $cds_end && @exon_start_arr>0){
			annotate_transc_exon_frames("$$id_name",	$contig, ($strand=~/-/?-1:1), $cds_start, $cds_end, \@exon_start_arr, \@exon_end_arr);
		}
		$cds_start="";
		$cds_end="";
	}	
}
close GENE_ANO;

##2.1, read junction_inf and load hash %gene2sample2juncNum and %gene2sample2juncReadNum (calculate average counts for junctions in a gene in a sample)
open (INJUNC,$junction_inf) || die "error open $junction_inf\n";
print "#open $junction_inf\n";
print LOGF "#open $junction_inf\n";
my @readNum_names=();
while(<INJUNC>){
	chomp;
	if(/juncpos5/){ #readType refseqid  row  juncpos5 juncpos3 num_DMSO_rep1  num_DMSO_rep2   num_DMSO_rep3   num_ro247_300nM_rep1 ... gene_symbol  gene_id contig  strand  gene_desc
		@headername_arr=split(/\t/);
		print $_."\n";
		if($sampleNames_str){
			@readNum_names=map("num_$_", split(/\s+/,$sampleNames_str) );
		}elsif($sampleNames_autoDetect_patt){
			foreach my $head_name(@headername_arr){
				if($head_name=~/$sampleNames_autoDetect_patt/){
					push @readNum_names, $head_name;
				}
			}
		}
		if(scalar @readNum_names>0){
			print "\@readNum_names=".join(" ",@readNum_names).";\n";
		}else{
			print "Warning: \@readNum_names not defined !\n";
		}

	}else{
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}

		#calculate read number, 
		my $if_junc_expressed=1;
		if(scalar @readNum_names>0){
			map($gene2sample2juncNum{"$contig	$strand	$$id_name	$row"}{$_}+=($$_>0?1:0), @readNum_names); #total detected junction number in a gene in a sample
			map($gene2sample2juncReadNum{"$contig	$strand	$$id_name	$row"}{$_}+=$$_, @readNum_names); #total junction counts in a gene in a sample
		}
		
	}
}
close INJUNC;



##2.2, read junction_inf and load novel ss information; add new ss to %geneCut_point_h; load %all_transcript_units_stats_h; load %junc2totalReadNum_h
open (INJUNC,$junction_inf) || die "error open $junction_inf\n";
print "#open $junction_inf\n";
print LOGF "#open $junction_inf\n";
my @readNum_names=();
my %junc2totalReadNum_h=();
while(<INJUNC>){
	chomp;
	if(/juncpos5/){ #readType refseqid  row  juncpos5 juncpos3 num_DMSO_rep1  num_DMSO_rep2 ... gene_symbol  gene_id contig  strand  gene_desc
		@headername_arr=split(/\t/);
		print $_."\n";
		if($sampleNames_str){
			@readNum_names=map("num_$_", split(/\s+/,$sampleNames_str) );
		}elsif($sampleNames_autoDetect_patt){
			foreach my $head_name(@headername_arr){
				if($head_name=~/$sampleNames_autoDetect_patt/){
					push @readNum_names, $head_name;
				}
			}
		}
		if(scalar @readNum_names>0){
			print "\@readNum_names=".join(" ",@readNum_names).";\n";
		}else{
			print "Warning: \@readNum_names not defined !\n";
		}

	}else{
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}
		$counts{"num.junction.1.detected"}++;

		$all_juncs{"$contig	$strand	$juncpos5	$juncpos3"}++;
		

		#calculate read number, filter lowly detected novel junctions
		my $gene_info_id="$contig	$strand	$$id_name	$row";
		my $if_junc_expressed=1;
		if(scalar @readNum_names>0){
			my $total_read_num=sum(map($$_, @readNum_names));
			my $total_sample_num=sum(map($$_>0?1:0, @readNum_names));
			$junc2totalReadNum_h{"$contig	$strand	$juncpos5	$juncpos3"}=$total_read_num;
			# foreach my $readNum_name1(@readNum_names){
			# 	if($$readNum_name1>10) {
			# 		print "$gene_info_id: $readNum_name1, $$readNum_name1 > $gene2sample2juncReadNum{$gene_info_id}{$readNum_name1}/$gene2sample2juncNum{$gene_info_id}{$readNum_name1}*$novelJunc_SuppRelReadMin; \n";	
			# 	}
			# }; #debug
			my $SampleNum_GreaterThan_RelCounts = sum(map(
				$$_ > $gene2sample2juncReadNum{$gene_info_id}{$_}/ ($gene2sample2juncNum{$gene_info_id}{$_}? $gene2sample2juncNum{$gene_info_id}{$_}:1) * $novelJunc_SuppRelReadMin ? 1:0
				,@readNum_names ) );
			$if_junc_expressed=$total_read_num>=$novelJunc_SuppReadMin && $total_sample_num>=$novelJunc_SuppSampleMin && $SampleNum_GreaterThan_RelCounts>1 ;
		}
		$counts{"num.junction.2.expressed (total read number>=$novelJunc_SuppReadMin & samples detected>=$novelJunc_SuppSampleMin & Counts/Average>$novelJunc_SuppRelReadMin)"}++ if $if_junc_expressed;
		if ($if_junc_expressed){
			$all_transcript_units_stats_h{$gene_info_id}++ ;
			$gene2junc_h{$gene_info_id}{"$juncpos5	$juncpos3"}++;
		}
		

		my (%junc_ano_h);
		foreach my $ssType('5','3'){
			my $ss_key="$contig:$strand:".${"juncpos$ssType"};
			if($kn_ss_h{$ssType}{$ss_key}{$$id_name}){ #known ss
				$junc_ano_h{$ssType}=$kn_ss_h{$ssType}{$ss_key}{$$id_name};
			}else{
				#search if ss in other transcript ids
				if(scalar keys %{$kn_ss_h{$ssType}{$ss_key}}>1 ){
					foreach my $other_id(keys %{$kn_ss_h{$ssType}{$ss_key}}){
						next if $other_id eq $$id_name;
						$junc_ano_h{$ssType}="$other_id:".$kn_ss_h{$ssType}{$ss_key}{$other_id};
						last;
					}
				}else{
					$junc_ano_h{$ssType}="n";
				}
				if($junc_ano_h{$ssType} ne "n" || $if_junc_expressed){
					$geneCut_point_h{$gene_info_id}{${"juncpos$ssType"}}="ss$ssType	".$junc_ano_h{$ssType};
				}
			}

		} #end foreach my $ssType
		if ($if_junc_expressed){
			my $if_junc_novel = $junc_ano_h{'5'} eq "n" && $junc_ano_h{'3'} eq "n";
			$counts{"num.junction.3.".($if_junc_novel?"novel":"known")}++;
		}
		
	}
}
close INJUNC;

$counts{"num.gene.1.detected"}=scalar keys %all_transcript_units_stats_h;



#3, output all junction information to $outAllJunc_f
open (OUTALLJUNC,">$outAllJunc_f") || die "error write $outAllJunc_f\n";
print OUTALLJUNC "contig	strand	juncpos5	juncpos3\n";
foreach my $junc_info (sort keys %all_juncs){
	print OUTALLJUNC "$junc_info\n";
}
close OUTALLJUNC;


#4, for each gene/transcript, cut into pieces based on known and novel splice site, transcript start and end, CDS start and end; output
open (OUTFLAT_F,">$outGeneBlock_flatf" ) || die "error write $outGeneBlock_flatf\n";
open (OUTANO_F,">$outGeneBlock_anof" ) || die "error write $outGeneBlock_anof\n";
open (BLO_BOUND_ANO_F,">$outBlockBoundary_f" ) || die "error write $outBlockBoundary_f\n";
print OUTANO_F "contig	strand	$id_name	row	Gblock_id	start_pos	end_pos	startPosType	startPosAno	endPosType	endPosAno	reg_len	region_ano\n";
print BLO_BOUND_ANO_F "contig	strand	$id_name	row	site_pos	PosType	PosAno\n";
print OUTFLAT_F "contig	strand	$id_name	row	Gblock_id	transc_start	transc_end	exon_starts	exon_ends\n";
print "#write $outGeneBlock_flatf\n#write $outGeneBlock_anof\n";


my %as2outhandle_h=();
if($outAS_anofroot){
	foreach my $as_type(@AS_types){
		my $as_ano_out_f="$outAS_anofroot$as_type.txt";
		my $as_ano_out_fh="OUTFILE_HANDLE.$as_ano_out_f";
		open ($as_ano_out_fh, ">$as_ano_out_f") || die "error write $as_ano_out_f\n";
		print "#write $as_ano_out_f\n";
		print $as_ano_out_fh "contig	strand	$id_name	row	as_id	". $asType2_anoHeaders_h{$as_type}."\n";
		$as2outhandle_h{$as_type}=$as_ano_out_fh;
	}
}

foreach my $transcript_unit_str (sort keys %all_transcript_units_stats_h){
	my ($contig,$strand,$id,$row)=split(/\t/,$transcript_unit_str);
	my @cut_points=sort {$a<=>$b} keys %{$geneCut_point_h{$transcript_unit_str}};
	my $strand_s=1;
	my %onegene_junc5to3_h=();
	my %onegene_junc3to5_h=();
	my %onegene_exonReg_Start2End_h=();
	my %onegene_exonReg_End2Start_h=();

	foreach my $juncs (keys %{$gene2junc_h{$transcript_unit_str}} ){
		my ($ss5,$ss3)=split(/\t/, $juncs);
		#print "\$juncs=$juncs; $ss5,$ss3\n";
		$onegene_junc5to3_h{$ss5}{$ss3}=1;
		$onegene_junc3to5_h{$ss3}{$ss5}=1;
	}
	if($strand eq "-"){
		$strand_s=-1;
		@cut_points=reverse(@cut_points);
	}

	my $current_kn_exon_Ano="";
	foreach my $cut_pos(@cut_points){ #write BLO_BOUND_ANO_F
		my $pos_info=$geneCut_point_h{$transcript_unit_str}{$cut_pos};
		my ($PosType,$PosAno)=split(/\t/,$pos_info);
		if($PosAno eq "n" && $current_kn_exon_Ano){ $PosAno="n-d$current_kn_exon_Ano"; } #for novel ss, annotate as closest upstream exon id (n-d means novel-downstream ...)
		if($PosAno !~/^n-d/){
			$current_kn_exon_Ano="$PosAno$PosType";
		}
		print BLO_BOUND_ANO_F "$transcript_unit_str	$cut_pos	$PosType	$PosAno\n";
	}

	my $block_id=0;
	my $current_kn_exon_Ano="";
	for(my $i=0;$i<@cut_points-1;$i++){
		my $start_pos=$cut_points[$i];
		my $end_pos=$cut_points[$i+1];
		my $start_pos_info=$geneCut_point_h{$transcript_unit_str}{$start_pos};
		my $end_pos_info=$geneCut_point_h{$transcript_unit_str}{$end_pos};
		my ($startPosType,$startPosAno)=split(/\t/,$start_pos_info);
		my ($endPosType,$endPosAno)=split(/\t/,$end_pos_info);
		if($startPosAno eq "n" && $current_kn_exon_Ano){ $startPosAno="n-d$current_kn_exon_Ano"; } #for novel ss, annotate as closest upstream exon id (n-d means novel-downstream ...)
		if($endPosAno eq "n" && $current_kn_exon_Ano){ $endPosAno="n-d$current_kn_exon_Ano"; }
		if($endPosAno !~/^n-d/){
			$current_kn_exon_Ano="$endPosAno$endPosType";
		}elsif($startPosAno !~/^n-d/ ){
			$current_kn_exon_Ano="$startPosAno$startPosType";
		}


		#modify position
		if($startPosType=~/te|ce|ss5/ ){
			$start_pos+=$strand_s;
		}
		if($endPosType=~/ts|cs|ss3/ ){
			$end_pos-=$strand_s;
		}
		my $reg_len=abs($start_pos-$end_pos)+1;
		my $region_ano=annotate_exonic_region($startPosType, $endPosType, $reg_len);
		
		if($region_ano=~/exon/){
			$onegene_exonReg_Start2End_h{$start_pos}=$end_pos; #all these are exonic region 1-based coordinates
			$onegene_exonReg_End2Start_h{$end_pos}=$start_pos;
		}

		if($start_pos*$strand_s>$end_pos*$strand_s){ #error
			print LOGF "#error found: $transcript_unit_str; start_pos=$start_pos($start_pos_info); end_pos=$end_pos($end_pos_info)\n";
			$counts{"num.gene_subRegion.1.error"}++;
			next;
		}
		#annotate region
		
		#write output file:
		$block_id++;
		my $gblock_id="$id:$row:$block_id";
		print OUTANO_F "$transcript_unit_str	$gblock_id	$start_pos	$end_pos	$startPosType	$startPosAno	$endPosType	$endPosAno	$reg_len	$region_ano\n";
		my ($left_pos,$right_pos)=sort {$a<=>$b} ($start_pos, $end_pos);
		$left_pos-=1; #change to 0-based
		print OUTFLAT_F "$transcript_unit_str	$gblock_id	$left_pos	$right_pos	$left_pos	$right_pos\n"; 
	}

	#for a gene, find alternative splicing events:
	if($outAS_anofroot){
		foreach my $as_type(@AS_types){
			$subroutine_name="find_$as_type";
			my (@events)=&$subroutine_name($id, $contig, $strand_s, \%onegene_junc5to3_h, \%onegene_junc3to5_h, \%onegene_exonReg_Start2End_h, \%onegene_exonReg_End2Start_h);
			if(@events>0){
				#print AS output file
				my $as_ano_out_fh=$as2outhandle_h{$as_type};
				#print "$transcript_unit_str\n";
				#print join("\n",@events)."\n";
				print $as_ano_out_fh join("", map("$transcript_unit_str	".($_+1)."	$events[$_]\n",0..(@events-1)) );
				$counts{"num.AS.out.$as_type"}+=scalar @events;
			}
		}
	}

}
close OUTFLAT_F;
close OUTANO_F;
close BLO_BOUND_ANO_F;

if($outAS_anofroot){
	foreach my $as_type(@AS_types){
		my $as_ano_out_fh=$as2outhandle_h{$as_type};
		close $as_ano_out_fh;
	}
}



print  "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
print LOGF "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
close LOGF;


sub annotate_exonic_region{
	my ($startPosType, $endPosType, $reg_len)=@_;
	my $region_ano="";
	if($startPosType=~/ss3/ && $endPosType=~/ss5/ ){
		$region_ano="exon";
	}elsif($startPosType=~/ss5/ && $endPosType=~/ss3/){
		$region_ano="intron";
	}elsif($startPosType=~/ss3/ && $endPosType=~/ss3/){
		if($reg_len<200){
			$region_ano="exon_a3ss";
		}else{
			$region_ano="intron_a3ss";
		}
	}elsif($startPosType=~/ss5/ && $endPosType=~/ss5/){
		if($reg_len<200){
			$region_ano="exon_a5ss";
		}else{
			$region_ano="intron_a5ss";
		}
	}elsif($startPosType=~/cs|ce/ || $endPosType=~/cs|ce/){
		$region_ano="exon_part";
	}elsif($startPosType=~/ss3/ || $endPosType=~/ss5/){
		$region_ano="exon_part";
	}elsif($startPosType=~/ss5/ || $endPosType=~/ss3/){
		$region_ano="intron_part";
	}

	$region_ano;
}

sub find_a5ss(){
	my($transcript_id, $contig, $strand_s, $onegene_junc5to3_h_addr, $onegene_junc3to5_h_addr, $onegene_exonReg_Start2End_h_addr, $onegene_exonReg_End2Start_h_addr)=@_;
	my %onegene_junc5to3_h=%$onegene_junc5to3_h_addr;
	my %onegene_junc3to5_h=%$onegene_junc3to5_h_addr;
	my %onegene_exonReg_Start2End_h=%$onegene_exonReg_Start2End_h_addr;
	my %onegene_exonReg_End2Start_h=%$onegene_exonReg_End2Start_h_addr;
	my $strand=$strand_s == -1? "-":"+";
	#print "find_a5ss $strand_s\n";
	my @result_arr=();
	my @aa3s=sort {$a<=>$b} keys %onegene_junc3to5_h;
	if($strand_s == -1){ @aa3s=reverse(@aa3s); }
	foreach my $ss3(@aa3s){
		my @ss5s=sort {$a<=>$b} keys %{$onegene_junc3to5_h{$ss3}};
		if($strand_s == -1){ @ss5s=reverse(@ss5s); }
		#print "ss3=$ss3; ss5=". join(" ",@ss5s). ";\n";
		next if scalar @ss5s <=1;
		for( my $i=0;$i<=(@ss5s-2);$i++){
			if( $onegene_exonReg_Start2End_h{ $ss5s[$i]+$strand_s } eq $ss5s[$i+1]  ){ #find alternative exonic region 
				my $exon_i=($ss5s[$i]+$strand_s)."-".$ss5s[$i+1] ; #alternative exonic region
				my $exon_uc=$onegene_exonReg_End2Start_h{$ss5s[$i]}."-".$ss5s[$i]  ; #common region of upstream exon
				my $exon_dc= $ss3."-".$onegene_exonReg_Start2End_h{$ss3}; ##common region of downstream exon
				my $junc_e=$ss5s[$i]."-".$ss3; #junction resulting for exclusion of alternative region
				my $junc_i=$ss5s[$i+1]."-".$ss3; #junction resulting for inclusion of alternative region
				my $allJunc_ReadNum=$junc2totalReadNum_h{"$contig	$strand	".$ss5s[$i]."	$ss3"} + $junc2totalReadNum_h{"$contig	$strand	".$ss5s[$i+1]."	$ss3"};

				#annotate frame information for exon_i if possible
				my $exon_i_len=abs(($ss5s[$i]+$strand_s)-$ss5s[$i+1]) +1;
				my $exon_i_inFrame="NA";
				
				my $exon_i_startFrame="NA";
				my $exon_i_seq=uc( extract_seq_fr_contig( $contig_file_dir, $contig, $strand_s, $exon_i, $contig_f  ) );
				my $exon_i_InFrameSeq=$exon_i_seq;
				my $exon_i_aa=""; my $exon_i_stopCodon="NA";
				my ($frame,$part_seq,$cum_CDS_length);
				my $use_transcript_id="";
				if($ss5_frame_ano_h{"$contig:$strand_s:".$ss5s[$i]}{$transcript_id} || $ss5_frame_ano_h{"$contig:$strand_s:".$ss5s[$i+1]}{$transcript_id} ){ #   || $ss3_frame_ano_h{ "$contig:$strand_s:".$ss3 }{$transcript_id}
					$use_transcript_id=$transcript_id;  #the representative transcript has coding annotation
				}else{ #annotation may come from other transcript in the same gene
					my @use_transcript_ids= keys %{ $ss5_frame_ano_h{"$contig:$strand_s:".$ss5s[$i]} };
					if(scalar @use_transcript_ids == 0){
						@use_transcript_ids= keys %{ $ss5_frame_ano_h{"$contig:$strand_s:".$ss5s[$i+1]} };
						#if(scalar @use_transcript_ids == 0){
						#	@use_transcript_ids= keys %{ $ss3_frame_ano_h{ "$contig:$strand_s:".$ss3 } };
						#}
					}
					if(@use_transcript_ids>0){
						$use_transcript_id=shift @use_transcript_ids;
					}
				}
				if($use_transcript_id){ #if there are coding annotation for the gene
					($frame,$part_seq,$cum_CDS_length)=split(/\t/, $ss5_frame_ano_h{"$contig:$strand_s:".$ss5s[$i]}{$use_transcript_id} ) ;
					if( $frame=~/^[012]/ ){ #frame from upstream ss5
						$exon_i_InFrameSeq=lc($part_seq).$exon_i_InFrameSeq;
						$exon_i_seq=lc($part_seq).$exon_i_seq;
					}elsif($exon_i_inFrame eq "Y"){ #frame from downstream 5'ss or 3'ss, annotate only when the alternative region is in-frame
						($frame,$part_seq,$cum_CDS_length) = split(/\t/, $ss5_frame_ano_h{"$contig:$strand_s:".$ss5s[$i+1]}{$use_transcript_id} );
						#if( $frame !~ /^[012]/ ){ 
						#	($frame,$part_seq,$cum_CDS_length) = split(/\t/, $ss3_frame_ano_h{ "$contig:$strand_s:".$ss3 }{$use_transcript_id} );
						#}
						if($frame=~/^[012]/ && $frame>0){
							$exon_i_InFrameSeq=substr($exon_i_InFrameSeq,(3-$frame),-$frame);
						}

					}
					if($frame ne "NA" && $frame ne ""){
						$exon_i_inFrame=$exon_i_len % 3 eq '0' ? "Y":"N";
						$exon_i_aa=translate_one_DNA($exon_i_InFrameSeq);
						$exon_i_stopCodon=$exon_i_aa=~/-/ ? "Y":"N";
					}
				}


				push @result_arr, "$exon_i	$exon_uc	$exon_dc	$junc_i	$junc_e	$exon_i_len	$exon_i_inFrame	$exon_i_seq	$exon_i_aa	$exon_i_stopCodon	$allJunc_ReadNum";
			}
		}
	}
	return (@result_arr);
}

sub find_a3ss(){
	my($transcript_id, $contig, $strand_s, $onegene_junc5to3_h_addr, $onegene_junc3to5_h_addr, $onegene_exonReg_Start2End_h_addr, $onegene_exonReg_End2Start_h_addr)=@_;
	my %onegene_junc5to3_h=%$onegene_junc5to3_h_addr;
	my %onegene_junc3to5_h=%$onegene_junc3to5_h_addr;
	my %onegene_exonReg_Start2End_h=%$onegene_exonReg_Start2End_h_addr;
	my %onegene_exonReg_End2Start_h=%$onegene_exonReg_End2Start_h_addr;
	my @result_arr=();
	my $strand=$strand_s == -1? "-":"+";
	my @aa5s=sort {$a<=>$b} keys %onegene_junc5to3_h;
	if($strand_s == -1){ @aa5s=reverse(@aa5s); }	
	foreach my $ss5(@aa5s){
		my @ss3s=sort {$a<=>$b} keys %{$onegene_junc5to3_h{$ss5}};
		if($strand_s == -1){ @ss3s=reverse(@ss3s); }
		next if scalar @ss3s <=1;
		for( my $i=0;$i<=(@ss3s-2);$i++){
			if( $onegene_exonReg_Start2End_h{ $ss3s[$i] } eq $ss3s[$i+1]-$strand_s  ){ #find alternative exonic region 
				my $exon_i=($ss3s[$i])."-".($ss3s[$i+1]-$strand_s) ; #alternative exonic region
				my $exon_uc=$onegene_exonReg_End2Start_h{$ss5}."-".$ss5  ; #common region of upstream exon
				my $exon_dc= $ss3s[$i+1]."-".$onegene_exonReg_Start2End_h{$ss3s[$i+1]}; ##common region of downstream exon
				my $junc_e=$ss5."-".$ss3s[$i+1]; #junction resulting for exclusion of alternative region
				my $junc_i=$ss5."-".$ss3s[$i]; #junction resulting for inclusion of alternative region
				my $allJunc_ReadNum=$junc2totalReadNum_h{"$contig	$strand	$ss5	".$ss3s[$i]} + $junc2totalReadNum_h{"$contig	$strand	$ss5	".$ss3s[$i+1]} ;
				#annotate frame information for exon_i if possible
				my $exon_i_len=abs($ss3s[$i]-($ss3s[$i+1]-$strand_s)) +1;
				my $exon_i_inFrame="NA";
				
				my $exon_i_startFrame="NA";
				my $exon_i_seq=uc( extract_seq_fr_contig( $contig_file_dir, $contig, $strand_s, $exon_i, $contig_f  ) );
				my $exon_i_InFrameSeq=$exon_i_seq;
				my $exon_i_aa=""; my $exon_i_stopCodon="NA";
				my ($frame,$part_seq,$cum_CDS_length);
				my $use_transcript_id="";
				if($ss5_frame_ano_h{"$contig:$strand_s:".$ss5}{$transcript_id} || $ss3_frame_ano_h{ "$contig:$strand_s:".$ss3s[$i] }{$transcript_id}){
					$use_transcript_id=$transcript_id;  #the representative transcript has coding annotation
				}else{ #annotation may come from other transcript in the same gene
					my @use_transcript_ids= keys %{$ss5_frame_ano_h{"$contig:$strand_s:".$ss5}};
					if(@use_transcript_ids>0){
						$use_transcript_id=shift @use_transcript_ids;
					}else{
						@use_transcript_ids= keys %{$ss3_frame_ano_h{ "$contig:$strand_s:".$ss3s[$i] }};
						if(@use_transcript_ids>0){  $use_transcript_id=shift @use_transcript_ids; }
					}
				}
				if($use_transcript_id){ #if there are coding annotation for the gene
					($frame,$part_seq,$cum_CDS_length)=split(/\t/, $ss5_frame_ano_h{"$contig:$strand_s:".$ss5}{$use_transcript_id} ) ;
					if( $frame=~/^[012]/ ){
						$exon_i_InFrameSeq=lc($part_seq).$exon_i_InFrameSeq;
						$exon_i_seq=lc($part_seq).$exon_i_seq;
					}else{
						($frame,$part_seq,$cum_CDS_length) = split(/\t/, $ss3_frame_ano_h{ "$contig:$strand_s:".$ss3s[$i] }{$use_transcript_id} );
						if( $frame=~/^[012]/ ){
							$exon_i_InFrameSeq=substr($exon_i_InFrameSeq, $frame eq 0? 0: (3-$frame) );
						}
					}
					if($frame ne "NA" && $frame ne ""){
						$exon_i_inFrame=$exon_i_len % 3 eq '0' ? "Y":"N";
						$exon_i_aa=translate_one_DNA($exon_i_InFrameSeq);
						$exon_i_stopCodon=$exon_i_aa=~/-/ ? "Y":"N";
					}
				}

				push @result_arr, "$exon_i	$exon_uc	$exon_dc	$junc_i	$junc_e	$exon_i_len	$exon_i_inFrame	$exon_i_seq	$exon_i_aa	$exon_i_stopCodon	$allJunc_ReadNum";
			}
		}
	}
	return (@result_arr);
}


sub find_se(){
	my($transcript_id, $contig, $strand_s, $onegene_junc5to3_h_addr, $onegene_junc3to5_h_addr, $onegene_exonReg_Start2End_h_addr, $onegene_exonReg_End2Start_h_addr)=@_;
	my %onegene_junc5to3_h=%$onegene_junc5to3_h_addr;
	my %onegene_junc3to5_h=%$onegene_junc3to5_h_addr;
	my %onegene_exonReg_Start2End_h=%$onegene_exonReg_Start2End_h_addr;
	my %onegene_exonReg_End2Start_h=%$onegene_exonReg_End2Start_h_addr;
	my @result_arr=();
	my @aa5s=sort {$a<=>$b} keys %onegene_junc5to3_h;
	if($strand_s == -1){ @aa5s=reverse(@aa5s); }
	my $strand=$strand_s == -1? "-":"+";
	foreach my $ss5(@aa5s){
		my @ss3s=sort {$a<=>$b} keys %{$onegene_junc5to3_h{$ss5}};
		if($strand_s == -1){ @ss3s=reverse(@ss3s); }
		next if scalar @ss3s <=1;
		for( my $i=0;$i<=(@ss3s-2);$i++){
			for(my $j=$i+1;$j<=(@ss3s-1);$j++){
				my $se_end=$onegene_exonReg_Start2End_h{ $ss3s[$i] };

				if( $onegene_junc5to3_h{$se_end}{$ss3s[$j]}  ){
					my $exon_i=($ss3s[$i])."-".$se_end ; #alternative exonic region
					my $exon_uc=$onegene_exonReg_End2Start_h{$ss5}."-".$ss5  ; #common region of upstream exon
					my $exon_dc= $ss3s[$j]."-".$onegene_exonReg_Start2End_h{$ss3s[$j]}; ##common region of downstream exon
					my $junc_e=$ss5."-".$ss3s[$j]; #junction resulting for exclusion of alternative region
					my $junc_ui=$ss5."-".$ss3s[$i]; #upstream junction resulting for inclusion of alternative region
					my $junc_di=$se_end."-".$ss3s[$j]; #downstream junction resulting for inclusion of alternative region
					
					my $allJunc_ReadNum=$junc2totalReadNum_h{"$contig	$strand	$ss5	".$ss3s[$j]} + $junc2totalReadNum_h{"$contig	$strand	$ss5	".$ss3s[$i]} + $junc2totalReadNum_h{"$contig	$strand	$se_end	".$ss3s[$j]};

					#annotate frame information for exon_i if possible
					my $exon_i_len=abs($ss3s[$i]-$se_end)+1;
					my $exon_i_inFrame="NA";

					my $exon_i_startFrame="NA";
					my $exon_i_seq=uc( extract_seq_fr_contig( $contig_file_dir, $contig, $strand_s, $exon_i, $contig_f  ) );
					my $exon_i_InFrameSeq=$exon_i_seq;
					my $exon_i_aa=""; my $exon_i_stopCodon="NA";
					my ($frame,$part_seq,$cum_CDS_length);
					my $use_transcript_id="";
					if($ss5_frame_ano_h{"$contig:$strand_s:".$ss5}{$transcript_id} || $ss3_frame_ano_h{ "$contig:$strand_s:".$ss3s[$i] }{$transcript_id}){
						$use_transcript_id=$transcript_id;  #the representative transcript has coding annotation
					}else{ #annotation may come from other transcript in the same gene
						my @use_transcript_ids= keys %{$ss5_frame_ano_h{"$contig:$strand_s:".$ss5}};
						if(@use_transcript_ids>0){
							$use_transcript_id=shift @use_transcript_ids;
						}else{
							@use_transcript_ids= keys %{$ss3_frame_ano_h{ "$contig:$strand_s:".$ss3s[$i] }};
							if(@use_transcript_ids>0){  $use_transcript_id=shift @use_transcript_ids; }
						}
					}
					if($use_transcript_id){ #if there are coding annotation for the gene
						($frame,$part_seq,$cum_CDS_length)=split(/\t/, $ss5_frame_ano_h{"$contig:$strand_s:".$ss5}{$use_transcript_id} ) ;
						if( $frame=~/^[012]/ ){
							$exon_i_InFrameSeq=lc($part_seq).$exon_i_InFrameSeq;
							$exon_i_seq=lc($part_seq).$exon_i_seq;
						}else{
							($frame,$part_seq,$cum_CDS_length) = split(/\t/, $ss3_frame_ano_h{ "$contig:$strand_s:".$ss3s[$i] }{$use_transcript_id} );
							if( $frame=~/^[012]/ ){
								$exon_i_InFrameSeq=substr($exon_i_InFrameSeq, $frame eq 0? 0: (3-$frame) );
							}
						}
						if($frame ne "NA" && $frame ne ""){
							$exon_i_aa=translate_one_DNA($exon_i_InFrameSeq);
							$exon_i_inFrame=$exon_i_len % 3 eq '0' ? "Y":"N";
							$exon_i_stopCodon=$exon_i_aa=~/-/ ? "Y":"N";
						}
					}
					push @result_arr, "$exon_i	$exon_uc	$exon_dc	$junc_ui	$junc_di	$junc_e	$exon_i_len	$exon_i_inFrame	$exon_i_seq	$exon_i_aa	$exon_i_stopCodon	$allJunc_ReadNum";
				}#end for ss3=$i
			} #end for ss3=$j
		}
	}
	return (@result_arr);
}



sub annotate_transc_exon_frames{ # load %ss5_frame_ano_h and %ss3_frame_ano_h
	my ($transcript_id,	$contig, $strand_s, $cds_start, $cds_end, $exon_start_arr_r, $exon_end_arr_r)=@_;
	my @exon_start_arr=@$exon_start_arr_r;
	my @exon_end_arr=@$exon_end_arr_r;
	my $cum_CDS_length=0;
	foreach my $exon_i(1..@exon_start_arr){
		my ($start_frame,$start_seq, $end_frame,$end_seq)=("NA","","NA","");
		my $start=$exon_start_arr[$exon_i-1];
		my $end=$exon_end_arr[$exon_i-1];
		if(!$start || !$end){
			die "error in annotate_transc_exon_frames, exon coordinates not defined (for: $transcript_id)!\n";
		}
		my $if_cdsStart_inExon=judge_1pos_rel_1reg($start,$end,$cds_start,$strand_s) eq '0' ;
		my $if_cdsEnd_inExon=judge_1pos_rel_1reg($start,$end,$cds_end,$strand_s) eq '0';
		my $if_exonStart_inCDS=judge_1pos_rel_1reg($cds_start,$cds_end,$start,$strand_s) eq '0';
		my $if_exonEnd_inCDS=judge_1pos_rel_1reg($cds_start,$cds_end,$end,$strand_s) eq '0';
		#print "cds:$cds_start-$cds_end; exon:$start-$end; if_cdsStart_inExon=$if_cdsStart_inExon; if_cdsEnd_inExon=$if_cdsEnd_inExon; if_exonStart_inCDS=$if_exonStart_inCDS; if_exonEnd_inCDS=$if_exonEnd_inCDS\n";
		if( $if_cdsStart_inExon && $if_cdsEnd_inExon ){ #CDS start and end in a exon
			#nothing to annotate
		}else{
			if($if_exonStart_inCDS){
				$start_frame=$cum_CDS_length % 3;
				if($start_frame>0){ #extract remaining nts in start of a exon to form codon with upstream exon
					$start_seq=extract_seq_fr_contig($contig_file_dir, $contig, $strand_s, $start."-". ($start+(3-$start_frame-1)*$strand_s), $contig_f  );
				}
				$ss3_frame_ano_h{"$contig:$strand_s:$start"}{$transcript_id}="$start_frame	$start_seq	$cum_CDS_length";
			}
			if($if_exonEnd_inCDS){
				$cum_CDS_length+= abs($end- ($if_cdsStart_inExon?$cds_start:$start) )+1;
				$end_frame=$cum_CDS_length % 3;
				if($end_frame>0){
					$end_seq=extract_seq_fr_contig( $contig_file_dir, $contig, $strand_s, $end. "-". ($end-($end_frame-1)*$strand_s), $contig_f  );
				}
				$ss5_frame_ano_h{"$contig:$strand_s:$end"}{$transcript_id}="$end_frame	$end_seq	$cum_CDS_length";
			}

		}

	}
}

