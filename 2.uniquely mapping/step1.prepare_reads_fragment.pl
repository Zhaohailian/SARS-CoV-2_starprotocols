#!/usr/bin/perl
use strict;
die "perl $0 input.sam\n" if(@ARGV != 1);
my $input_sam=shift;
my $prefix=$input_sam;
$prefix=~s/.+\///;
$prefix=~s/.sam$//;

my $index;
open(FQ,">$prefix.fq") || die;
open(IS,$input_sam) || die;
while(my $line=<IS>){
	chomp $line;
	if($line=~/^@/){
		next;
	}
	my @sub=split/\t/,$line;
	$sub[1]=$sub[1] > 255 ? $sub[1]-256 : $sub[1];
	my ($skip,$match)=get_start_and_len($sub[5]);

	my $read_seq;
	$read_seq=substr($sub[9],$skip,$match);
	my $qual_seq;
	$qual_seq=substr($sub[10],$skip,$match);

	$index++;
	print FQ "@",$index,"_",$sub[2],":",$sub[3],"-",$sub[3]+$match-1,"_",$sub[0],"\n";
	print FQ $read_seq,"\n";		
	print FQ "+\n";
	print FQ $qual_seq,"\n";


	#my $qual="F" x length($read_seq);
	#print FQ $qual,"\n";
	#exit;
}


sub get_start_and_len{
        my $raw_cigar=shift;
	my $block_start=0;
        my $only_block_len=0;
	my $cigar=$raw_cigar;
        while($cigar=~/(\d+)(\w)/g){
                my $num=$1;
                my $class=$2;
		#print $num,"\t",$class,"\taa\n";
                if($class eq "S" or $class eq "H"){
			$block_start+=$1;
                }
		else{
			last;
		}
	}
	my $cigar=$raw_cigar;
	while($cigar=~/(\d+)(\w)/g){
		my $num=$1;
		my $class=$2;
		#print $num,"\t",$class,"\taa\n";
		if($class eq "S" or $class eq "H"){
			next;
		}
                elsif($class eq "M" or $class eq "I"){
                        $only_block_len+=$1;
                }
                elsif($class eq "D"){
                        next;
                }
                else{
                        die;
                }
        }
	#warn $cigar,"\t",$block_start,"\t",$only_block_len,"\n";
        return ($block_start,$only_block_len);
}


	
