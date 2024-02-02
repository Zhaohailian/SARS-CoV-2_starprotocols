#!/usr/bin/perl
die "perl $0 HeLa_merge.Chimeric.fq HeLa_merge.Chimeric.sam HeLa_merge.Chimeric.reMapHisat2.sam reMap_bowtie2_sam reMap_star_sam\n" if(@ARGV != 5);
my $input_chimeric_fq=shift;
my $input_chimeric_sam=shift;
my $reMap_hisat2_sam=shift;
my $reMap_bowtie2_sam=shift;
my $reMap_star_sam=shift;
my $outprefix=$input_chimeric_sam;
$outprefix=~s/.Chimeric.sam//;

my %badReads_basedOn_bowtie2;
my %goodReads_basedOn_bowtie2;
open(RBS,$reMap_bowtie2_sam) || die;
while(my $line=<RBS>){
	chomp $line;
	if($line=~/^@/){
		next;
	}
	my @sub=split/\s+/,$line;
	my $read_index=(split/_/,$sub[0],3)[0];
	my @read_input_loci1=(split/:/,$sub[0],2);
	my $read_input_loci=(split/_/,$read_input_loci1[1],2)[0];

	my $read_input_chr;
        my $read_input_start;
        my $read_input_end;

        if($read_input_loci=~/(\d+)-(\d+)/){
                $read_input_chr=(split/_/,$read_input_loci1[0],2)[1];
                $read_input_start=$1;
                $read_input_end=$2;
        }
        else{
                die;
        }

	#my $read_name=(split/_/,$sub[0],3)[2];
	my $read_name=(split/_/,$read_input_loci1[1],2)[1];
	if($sub[4] >= 20){
		my $read_remap_chr=$sub[2];
                my $read_remap_start=$sub[3];
                my $read_remap_end=$sub[3]+length($sub[9])-1;
		if($read_input_chr eq $read_remap_chr and $read_input_start < $read_remap_end and $read_input_end > $read_remap_start){
			$goodReads_basedOn_bowtie2{$read_index}=$read_name;
		}
		else{	
			$badReads_basedOn_bowtie2{$read_index}=$read_name;
		}
	}
	else{
		$badReads_basedOn_bowtie2{$read_index}=$read_name;
	}
}

my %badReads_basedOn_hisat2;
my %goodReads_basedOn_hisat2;
open(RHS,$reMap_hisat2_sam) || die;
while(my $line=<RHS>){
	chomp $line;
        if($line=~/^@/){
                next;
        }
        my @sub=split/\s+/,$line;
        my $read_index=(split/_/,$sub[0],3)[0];
	my @read_input_loci1=(split/:/,$sub[0],2);
	my $read_input_loci=(split/_/,$read_input_loci1[1],2)[0];
        my $read_input_chr;
        my $read_input_start;
        my $read_input_end;
        if($read_input_loci=~/(\d+)-(\d+)/){
                $read_input_chr=(split/_/,$read_input_loci1[0],2)[1];
                $read_input_start=$1;
                $read_input_end=$2;
        }
        else{
                die;
        }

	#my $read_name=(split/_/,$sub[0],3)[2];
	my $read_name=(split/_/,$read_input_loci1[1],2)[1];
	if($sub[4] >= 20){
		my $read_remap_chr=$sub[2];
                my $read_remap_start=$sub[3];
                my $read_remap_end=$sub[3]+length($sub[9])-1;
		if($read_input_chr eq $read_remap_chr and $read_input_start < $read_remap_end and $read_input_end > $read_remap_start){
			$goodReads_basedOn_hisat2{$read_index}=$read_name;
		}
		else{
			$badReads_basedOn_hisat2{$read_index}=$read_name;
		}
	}
	else{
		$badReads_basedOn_hisat2{$read_index}=$read_name;
     	}
}

my %badReads_basedOn_star;
my %goodReads_basedOn_star;
open(RSS,$reMap_star_sam) || die;
while(my $line=<RSS>){
        chomp $line;
        if($line=~/^@/){
                next;
        }
        my @sub=split/\s+/,$line;
        my $read_index=(split/_/,$sub[0],3)[0];
	my @read_input_loci1=(split/:/,$sub[0],2);
	my $read_input_loci=(split/_/,$read_input_loci1[1],2)[0];
        my $read_input_chr;
        my $read_input_start;
        my $read_input_end;
        if($read_input_loci=~/(\d+)-(\d+)/){
                $read_input_chr=(split/_/,$read_input_loci1[0],2)[1];
                $read_input_start=$1;
                $read_input_end=$2;
        }
        else{
                die;
        }
        #my $read_name=(split/_/,$sub[0],3)[2];
	my $read_name=(split/_/,$read_input_loci1[1],2)[1];
        if($sub[4] >= 20){
		my $read_remap_chr=$sub[2];
                my $read_remap_start=$sub[3];
                my $read_remap_end=$sub[3]+length($sub[9])-1;
		if($read_input_chr eq $read_remap_chr and $read_input_start < $read_remap_end and $read_input_end > $read_remap_start){
	                $goodReads_basedOn_star{$read_index}=$read_name;
		}
		else{
			$badReads_basedOn_star{$read_index}=$read_name;
		}
        }
        else{
                $badReads_basedOn_star{$read_index}=$read_name;
        }
}

####test####
#$goodReads_basedOn_bowtie2{16670182}=1;
#$goodReads_basedOn_hisat2{16670182}=1;
###test over

my $index_a=-1;
my $index_b=0;
my %firstHalf_uniq;
my %secondHalf_uniq;
open(OUTU,">$outprefix.bothUniq.1round.sam") || die;
open(OUTO,">$outprefix.Other.1round.sam") || die;
open(IS,$input_chimeric_sam) || die;
while(my $frag_a=<IS>){
	my $frag_b=<IS>;
	$index_a+=2;
	$index_b+=2;
	chomp $frag_a;
	chomp $frag_b;
	my @sub_a=split/\s+/,$frag_a;
	my @sub_b=split/\s+/,$frag_b;
	if((exists $goodReads_basedOn_bowtie2{$index_a} or exists $goodReads_basedOn_hisat2{$index_a} or exists $goodReads_basedOn_star{$index_a}) and (exists $goodReads_basedOn_bowtie2{$index_b} or exists $goodReads_basedOn_hisat2{$index_b} or exists $goodReads_basedOn_star{$index_b})){	#both end is uniq; by either method
		print OUTU $frag_a,"\n",$frag_b,"\n";
	}
	else{
		print OUTO $frag_a,"\n",$frag_b,"\n";
		if(exists $goodReads_basedOn_bowtie2{$index_a} and exists $goodReads_basedOn_hisat2{$index_a} and exists $goodReads_basedOn_star{$index_a}){	#first arm is uniq
			$firstHalf_uniq{$index_a}=1;
		}
		elsif(exists $goodReads_basedOn_bowtie2{$index_b} and exists $goodReads_basedOn_hisat2{$index_b} and exists $goodReads_basedOn_star{$index_b}){	#second arm is uniq
			$secondHalf_uniq{$index_b}=1;
		}
		else{
			next;
		}
	}
}

open(ICFQ,$input_chimeric_fq) || die;
open(UNIQ,">$outprefix.Chimeric.UniqOneArm.fq") || die;
open(MULT,">$outprefix.Chimeric.MultipleOneArm.fq") || die;
while(my $id_a=<ICFQ>){
	my $seq_a=<ICFQ>;
	my $symbol_a=<ICFQ>;
	my $qual_a=<ICFQ>;
	my $id_b=<ICFQ>;
	my $seq_b=<ICFQ>;
	my $symbol_b=<ICFQ>;
	my $qual_b=<ICFQ>;

	chomp $id_a;
	chomp $id_b;
	my $index_a=(split/_/,$id_a)[0];
	my $index_b=(split/_/,$id_b)[0];
	$index_a=~s/^@//;
	$index_b=~s/^@//;

	if(exists $firstHalf_uniq{$index_a} and exists $secondHalf_uniq{$index_b}){	#could not happen; should be classified as both uniq
		die;
	}
	elsif(exists $firstHalf_uniq{$index_a} and !exists $secondHalf_uniq{$index_b}){
		print UNIQ $id_a,"\n",$seq_a,$symbol_a,$qual_a;
		print MULT $id_b,"\n",$seq_b,$symbol_b,$qual_b;
	}
	elsif(!exists $firstHalf_uniq{$index_a} and exists $secondHalf_uniq{$index_b}){
		print MULT $id_a,"\n",$seq_a,$symbol_a,$qual_a;
		print UNIQ $id_b,"\n",$seq_b,$symbol_b,$qual_b;
	}
	else{	#both multi
		next;
	}
}








