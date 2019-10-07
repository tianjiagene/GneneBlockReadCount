#perl functions
#require "/Home/cocochen/analyze/perl.fun.inc.pl";
use POSIX;


sub create_dir_ifNotExist{ #input xxxx/xxxxx/  or xxx/xxx/xxx.xxx or /xxx/xxx/xxx
	$newdir=shift;
	$newdir=~s/[^\/]*$//;
	print "newdir=".$newdir."\n";
	if(!(-d $newdir)){
		@folder_arr=split(/\//,$newdir);
		my $i;
		for($i=0;$i<@folder_arr;$i++){
			$dir_i=join("/",@folder_arr[0..$i]);
			next if (!$dir_i); ##root
			if(! (-d $dir_i)){
				print "mkdir $dir_i, ";
				print(mkdir($dir_i,0774)."\n");
			}
		}
	}
}


sub extract_seq_fr_contig{ 
	my ($contig_file_dir,$contig_name,$strand,$regions_str,$contig_file)=@_;
	#if AACCCTTT + 3:6, return CCCT
	#if AACCCTTT - 3:6, return AGGG
	#$contig_name is the contig name following ">" in the contig file
	#when $search_reverse eq "T" and the strand is minus, program will search the reversed strand. the output will be reversed complementary NTs
	
	if(!$contig_seq{$contig_name} && !(scalar keys %contig_seq>0 && $contig_file) ){ #not loaded the seq
		#if (keys %contig_seq >6){%contig_seq={};}; #to clear memory, use this if the contig is very large
		$contig_file="$contig_file_dir/$contig_name.fa" if !$contig_file; 
		open (CONTIG_FILE_IN, $contig_file) || die ("error $contig_file\n"); #read the contig file and get the seq
		print "$contig_file opened...\n";
		my $contig_name_infile='';
		while (<CONTIG_FILE_IN>){
			s/\s+$//;
			if(substr($_,0,1)eq">"){
				$_=~s/>//;
				$contig_name_infile=$_;
				print " find sequence $contig_name_infile;\n";
				next;
			}
			$contig_seq{$contig_name_infile}.=$_;
		}
		close CONTIG_FILE_IN;
	}
	
	#extract seq from $contig_seq{$contig_name}
	my @regions=split(/,|\|/,$regions_str); #for multi-region (like junction)
	my $reg_id=0,$out_seq="";
	if($contig_seq{$contig_name}){
		while ($reg_id<@regions){
			my ($from,$to)=sort {$a<=>$b} split(/:|\-/,$regions[$reg_id]);
			my $seq_length=$to-$from+1;
			$out_seq.=substr($contig_seq{$contig_name},$from-1,$seq_length);
			$reg_id++;
		}
		if($strand=~/\-/){ # to reverse and complement
			$out_seq=reverse_complement($out_seq,"F");
		}
	}
	return uc($out_seq);
}

sub reverse_complement{ #for DNA or RNA sequence only
	my ($input_seq,$is_RNA)=@_;
	$input_seq=~tr/actguACTGU/tgacaTGACA/;
	if($is_RNA eq "T"){$input_seq=~tr/tT/uU/;}
	my $out_seq=reverse_str($input_seq);
	#for($i=length($input_seq)-1;$i>=0;$i--){
	#	$out_seq.=substr($input_seq,$i,1);
	#}
	return $out_seq;
}

sub reverse_str{
	my $input=shift;
	join("", reverse(split(//,$input)) )
}

sub return_time{
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	return (($year+1900)."/".($mon+1)."/".$mday." ".$hour.":".$min.":".$sec);
}



sub judge_2region_relation{ 
	#given 2 regions a and b (4 coordinates), output position of b relative to a (return =, 0, 1, -1, inc, ov1, ov-1)
	#if $ifout_overlapLen=1, output overlap info (overlap vode, overlap length, overlap coordinates and unoverlapped coordinates)
	my ($a1,$a2,$b1,$b2,$strand_a,$ifout_overlapLen)=@_;
	$strand_a=1 if ($strand_a eq "+");
	$strand_a=-1 if ($strand_a eq "-");
	$strand_a=1 if $strand_a ne '-1'; #$strand_a= 1 or -1 (represent + or - strand)
	if($a1>$a2){($a1,$a2)=($a2,$a1);  }
	($b1,$b2)=sort {$a<=>$b} ($b1,$b2);
	if(($a1 eq $b1 && $a2 eq $b2) || ($a1 eq $b2 && $a2 eq $b1) ){
		return $ifout_overlapLen ? ("=", $a2-$a1+1, "$a1:$a2", "") : "=";
	}else{
		$b1pos=judge_1pos_rel_1reg($a1,$a2,$b1);
		$b2pos=judge_1pos_rel_1reg($a1,$a2,$b2);
		if($b1pos==$b2pos){
			if($ifout_overlapLen){
				if($b1pos eq "0"){
					my @non_overlapped_a_arr=();
					if($b1>$a1){push @non_overlapped_a_arr, "$a1:".($b1-1);}
					if($b2<$a2){push @non_overlapped_a_arr, ($b2+1).":$a2";}
					return ($b1pos, abs($b2-$b1)+1,"$b1:$b2", join("|",@non_overlapped_a_arr) );
				}else{
					return ($b1pos*$strand_a, 0, "", "$a1:$a2");
				}
				
			}else{
				return $b1pos*$strand_a; #0 1 or -1
			}
			
		}else{
			if($b1pos*$b2pos==-1 || ($b1pos==-1 && $b2==$a2) ||  ($b2pos==1 && $b1==$a1) ){
				return $ifout_overlapLen ? ("inc", $a2-$a1+1, "$a1:$a2", "") : "inc"; #b include a
			}else{ #overlap
				my $posStr_rel=$b1pos+$b2pos; #if positive strand, relation: 1 or -1
				if($ifout_overlapLen){
					my @arr4=sort {$a<=>$b}($a1,$a2,$b1,$b2);
					my $overlap_len=$arr4[2]-$arr4[1]+1;
					return ("ov".($posStr_rel*$strand_a) , $overlap_len,  $arr4[1].":".$arr4[2],  $posStr_rel==1? "$a1:".($b1-1) : ($b2+1).":$a2"   );
				}else{
					return ( "ov".($posStr_rel*$strand_a) ); #ov1 or ov-1
				}
			}
		}		
	}
}


%match_syb2val_hash=(
	"="=>5,
	"0"=>4,
	"inc"=>3,
	"ov1"=>2,
	"ov-1"=>2,
	"1"=>1,
	"-1"=>1,
);

sub judge_1pos_rel_1reg{ #given 1 region a (a1 to a2) and another position x, judge position x relative to a (return -1, 0 or 1)
	my ($a1,$a2,$x,$strand_a)=@_; #a1 should<a2
	($a1,$a2)=sort{$a<=>$b}($a1,$a2);
	$strand_a=1 if ($strand_a eq "+");
	$strand_a=-1 if ($strand_a eq "-");
	$strand_a=1 if !$strand_a;
	if( ($a1-$x)*($a2-$x)<=0 ){
		return 0; #inside
	}else{
		if($x<$a1){
			return -1*$strand_a; #upstream if strand = 1 (+)
		}else{
			return $strand_a; #downstream if strand = 1 (+)
		}
	}
}

sub create_chr_reg_info_hash{ #create standard indexed chromosome region info hash
	my ($pos1,$pos2,$key1,$val_str,$chr_block_size,$key2)=@_; 
	#$key1 usally contains chromosome, [strand]
	#$key2 usally is a region (it should be unique for all values)
	my $posfrom=($pos1<$pos2)?$pos1:$pos2;
	my $posto=($pos1<$pos2)?$pos2:$pos1;
	my $posfrom_block_num=int($posfrom/$chr_block_size);
	my $posto_block_num=int($posto/$chr_block_size);
	$key2="$posfrom:$posto" if (!$key2);
	foreach my $pos_block_num(($posfrom_block_num..$posto_block_num)){
		$chro2reginfo_hash{"$key1:$pos_block_num"}{$key2}=$val_str;
	}
	($key1,$posfrom_block_num,$posto_block_num,$key2);
}


###############.sam parsing functions

sub cal_coverage_fr_cigar{ #calculate reads coverage region and junction/intron info from cigar string
	my ($from_pos,$ciger_str)=@_;
	return "" if $ciger_str!~/M/i;
	my @exonR_arr=(); #exon region
	my @intron_arr=();
	my $match_len;
	$ciger_str=~s/\d+[^\dMIDNSHP]//g;
	while($ciger_str){
		if($ciger_str=~s/^(\d+)M//){ #aligned region
			$match_len=$1;
			if(scalar @exonR_arr>0 && $exonR_arr[@exonR_arr-1] ne "|"){
				$exonR_arr[@exonR_arr-1]+=$match_len-1;
			}else{
				push (@exonR_arr, "$from_pos:", ($from_pos+$match_len-1) );
			}
			$from_pos+=$match_len-1;
		}
		if($ciger_str=~s/^(\d+)N//){ #gap, intron
			$match_len=$1;
			push (@intron_arr, ($exonR_arr[@exonR_arr-1]+1).":". ($exonR_arr[@exonR_arr-1]+$match_len) );
			push (@exonR_arr, "|") if @exonR_arr>0;
			$from_pos+=$match_len+1;
		}
		if($ciger_str=~s/^(\d+)[SHPI]//){
			#do nothing
		}
		if($ciger_str=~s/^(\d+)D//){ #deletion in reads
			$match_len=$1;
			$exonR_arr[@exonR_arr-1]+=$match_len+1 if @exonR_arr>0;
			$from_pos+=$match_len+1;
		}
	}
	return (join("",@exonR_arr), join("|",@intron_arr));
}


sub from_flag2strand { ###now only for tophat single end outfile
	my $flag_str=shift;
	($flag_str & 16)?"-":"+";
}

sub from_flag2info{
	my $flag_str=shift;
	#print "$flag_str	";
	my $if_PE=($flag_str & 1)?"PE":"SE";
	my $ifProperPair=($flag_str & 2)?1:0;
	my $ifUnmapped=($flag_str & 4)?1:0;
	my $Rstrand=($flag_str & 16)?"-":"+";
	my $mRstrand=$if_PE eq "PE"? (($flag_str & 0x0020)?"-":"+") : "";
	my $ifNotPrimary=($flag_str & 0x0100)?1:0;
	#print "$if_PE	$ifProperPair	$ifUnmapped	$Rstrand	$mRstrand	$ifNotPrimary\n";
	return($if_PE,$ifProperPair,$ifUnmapped,$Rstrand,$mRstrand,$ifNotPrimary);
}

sub judge_reads_quality{
	my ($OutSamFrom,$mismatch_cut,$clip_cut,$cigar,@tags_arr)=@_;
	##1, if clip length > clip_cut (bad quality)
	if($clip_cut ne 'none'){ #eg. set $clip_cut ='none' to invalid  this parameter
		my $clip_len_sum=0;
		while($cigar=~/(\d+)S/g){ # 1S22M
			$clip_len_sum+=$1;
			if($clip_len_sum>$clip_cut){
				return 0;
			}
		}
	}
	##2, if mismatch > $mismatch_cut (bad quality)
	if($mismatch_cut ne 'none'){
		my $tag_str=join("	",@tags_arr);
		my $mismatch=0;
		if($OutSamFrom eq "bwasw"){ #$cigar=1S22M  AS:i:22
			$cigar=~/(\d+)M/;
			my $match_len=$1;
			if($tag_str=~/AS:i:(\d+)/){
				$mismatch=($match_len-$1)/4 ; #1 mismatch score =-3, match score =1 for bwasw
			}
		}elsif($OutSamFrom=~/tophat|bowtie|STAR/){
			if($tag_str=~/NM:i:(\d+)/i){$mismatch=$1;} #NM:i:1
		}
		if($mismatch>$mismatch_cut){return 0;}else{return 1;}
	}
	return 1;
}

sub cigar2cover_len{
	my $cigar2=shift;
	$cigar2=~s/\d+[HSIP]//g;
	my $total_cover_len=0; #length in reference covered by reads aligned region
	foreach $len1(split(/[MND]/,$cigar2)){	$total_cover_len+=$len1;} #
	$total_cover_len;
}

sub cigar2read_len{
	my $cigar2=shift;
	$cigar2=~s/\d+[NDP]//g;
	my $readLen=0; #length in reference covered by reads aligned region
	foreach $len1(split(/[HSMI]/,$cigar2)){	$readLen+=$len1;} #
	$readLen;
}

###################
sub sum{
	my @invars=@_;
	my $sumval=0;
	map($sumval+=$_,@invars);
	$sumval;
}
sub average{
	my @invars=@_;
	return sum(@invars)/(scalar @invars) ;
}
sub stdev{
	my(@data) = @_;
	if(@data == 1){
	        return 0;
	}
	my $average = &average(@data);
	my $sqtotal = 0;
	foreach(@data) {
	        $sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@data-1)) ** 0.5;
	return $std;
}

sub sem{ #standard error of mean
	my(@data) = @_;
	return stdev(@data)/ (scalar @data ** 0.5) ;
}



%readsTypes2sco_h=( #higher score represent higher priority (applicable to RNA-seq); + - represent strand relative to reference gene only useful for directional reads 
	"CDS+"=>18,
	"CDS"=>17,
	"3UTR+"=>16,
	"3UTR"=>15,
	"5UTR+"=>14,
	"5UTR"=>13,
	"UTR3e+"=>12,
	"UTR3e"=>11,
	"intron+"=>10,
	"intron"=>9,
	"intron-"=>8,
	"UTR5e-"=>7,
	"UTR5e"=>6,
	"UTR5e+"=>5,
	"UTR3e-"=>4,
	"3UTR-"=>3,
	"CDS-"=>2,
	"5UTR-"=>1,
	"intergenic"=>0,
);

sub get_genomic_region{
	my ($ref_pos, $direction, $rel_pos1, $rel_pos2)=@_; #$direction = 1 or -1
	my $pos1=$ref_pos+$direction*$rel_pos1;
	my $pos2=$ref_pos+$direction*$rel_pos2;
	$pos1=1 if $pos1<=0;
	$pos2=1 if $pos2<=0;
	if($pos1>1 || $pos2>1){
		join(":", sort {$a<=>$b}($pos1,$pos2));
	}else{
		"";
	}
}


sub convert_regName2_coordinates{ #like: p101p2100, m100m60, convert to (101,2100) (-100,-60)
	my $region_name=shift; 
	$region_name=~s/([mp])/	$1/g;
	$region_name=~s/m/\-/g;
	$region_name=~s/p//g;
	$region_name=~s/^	//;
	split(/	/,$region_name);
}

sub bed2bigbed{
	my ($geno,$bedf,$group_name,$bed_name)=@_;
	my $bbf_name="$bedf.bb";
 	print "\n####convert $bedf to bigbed file: $bbf_name\n";
	system("sed \'1d\' $bedf | sort -k1,1 -k2,2n  > $bedf.sorted");
	system("bedToBigBed $bedf.sorted /Home/cocochen/soft/biosoft/genome/UCSC/$geno.chrom.sizes $bbf_name");
	system("rm $bedf.sorted");
	$bbf_name=~s/^.*\///;
	print  "####add track to UCSC:####\ntrack type=bigBed itemRgb=On visibility=3 colorByStrand=\'255,0,0 0,0,255\' group=$group_name priority=1 name=\"$bed_name\" description=\"$bed_name\" bigDataUrl=http://rna:rna\@rna.umdnj.edu/cocochen/bigwig/$group_name/$bbf_name\n\n";
}

sub translate_one_DNA {
	my $seq=@_[0];
	$seq=~s/\w{3}/$& /g;
	$seq=~s/ \w{1,2}$//g;
	$seq=~s/TTT|TTC/F/gi;
	$seq=~s/TAT|TAC/Y/gi;
	$seq=~s/CAT|CAC/H/gi;
	$seq=~s/CAA|CAG/Q/gi;
	$seq=~s/AAT|AAC/N/gi;
	$seq=~s/AAA|AAG/K/gi;
	$seq=~s/GAT|GAC/D/gi;
	$seq=~s/GAA|GAG/E/gi;
	$seq=~s/TGT|TGC/C/gi;
	$seq=~s/TTA|TTG|CTT|CTC|CTA|CTG/L/gi;
	$seq=~s/AGT|AGC|TCT|TCC|TCA|TCG/S/gi;
	$seq=~s/AGA|AGG|CGT|CGC|CGA|CGG/R/gi;
	$seq=~s/ATT|ATC|ATA/I/gi;
	$seq=~s/TGG/W/gi;
	$seq=~s/ATG/M/gi;
	$seq=~s/GTT|GTC|GTA|GTG/V/gi;
	$seq=~s/CCT|CCA|CCC|CCG/P/gi;
	$seq=~s/ACT|ACA|ACC|ACG/T/gi;
	$seq=~s/GCT|GCA|GCC|GCG/A/gi;
	$seq=~s/GGT|GGA|GGC|GGG/G/gi;
	$seq=~s/TAA|TAG|TGA/-/gi; # stop codon to -
	$seq=~s/\w{3}/?/g; #other unconverted
	$seq=~s/ //g;
	return $seq;
}

@all_codons=split(/;/,"TTT;TTC;TAT;TAC;CAT;CAC;CAA;CAG;AAT;AAC;AAA;AAG;GAT;GAC;GAA;GAG;TGT;TGC;TTA;TTG;CTT;CTC;CTA;CTG;AGT;AGC;TCT;TCC;TCA;TCG;AGA;AGG;CGT;CGC;CGA;CGG;ATT;ATC;ATA;TGG;ATG;GTT;GTC;GTA;GTG;CCT;CCA;CCC;CCG;ACT;ACA;ACC;ACG;GCT;GCA;GCC;GCG;GGT;GGA;GGC;GGG;TAA;TAG;TGA;");
sub count_codon_freq{
	my $seq=@_[0];
	$seq=~s/\w{3}/$& /g;
	$seq=uc($seq);
	my @codons=split(/ /, $seq);
	my %codon_fre=();
	map($codon_fre{$_}++, @codons);
	return( map(($codon_fre{$_}?$codon_fre{$_}:0),  @all_codons) );
}


sub cal_refseq_genelen{ #input are from ucsc genome browser refflat file (or other similar format)
	my ($cds_start,$cds_end,$exon_starts,$exon_ends)=@_;
	my @exon_start_arr=split(/,/,$exon_starts);
	my @exon_end_arr=split(/,/,$exon_ends);	
	my $cds_len=0; my $transcript_len=0;
	my $if_cds=0;
	for(my $i=0;$i<@exon_start_arr;$i++){ 
		my $exon_start=$exon_start_arr[$i];
		my $exon_end=$exon_end_arr[$i];
		$transcript_len+= abs($exon_end-$exon_start);
		
		my $cds_start2exon= judge_1pos_rel_1reg($exon_start,$exon_end, $cds_start);
		my $cds_end2exon= judge_1pos_rel_1reg($exon_start,$exon_end, $cds_end);
		if($cds_start2exon eq "0"){
			$exon_start=$cds_start;
			$if_cds=1;
		}
		if($cds_end2exon eq "0"){
			$exon_end=$cds_end;
			$cds_len += abs($exon_end-$exon_start);
			$if_cds=0;
		}
		
		if($if_cds){
			$cds_len+= abs($exon_end-$exon_start);
		}
	}	
	return ($cds_len,$transcript_len);
}

1;


