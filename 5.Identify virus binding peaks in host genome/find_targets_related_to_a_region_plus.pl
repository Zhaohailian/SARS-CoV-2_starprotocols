#!/usr/bin/perl

die "perl $0 <all_pairs.sam> <region_chr> <region_start> <region_end>\n" if(@ARGV != 4);
my $contact_sam=shift;
my $chr=shift;
my $start=shift;
my $end=shift;

open(CS,$contact_sam) || die;
while(my $frag_a=<CS>){
	if($frag_a=~/^@/){
		print $frag_a;
		next;
	}
	my $frag_b=<CS>;
	my @sub_a=split/\s+/,$frag_a;
	my @sub_b=split/\s+/,$frag_b;
	my $judge_a=good($frag_a);
	my $judge_b=good($frag_b);
	if($judge_a=~/good/ and $judge_b=~/good/){	#both
		next;
	}
	elsif($judge_a=~/good/ and $judge_b=~/bad/){
		print $frag_a;
		print $frag_b;
	}
	elsif($judge_a=~/bad/ and $judge_b=~/good/){
		print $frag_a;
		print $frag_b;
	}
	else{
		next;
	}
}
	


sub good{
	my $frag=shift;
	my @arr=split/\s+/,$frag;
	my @info=split/_/,$arr[0];
	if(($info[2] ne "Plus") and ($info[3] ne "Plus")){
		return "bad";
	}
	if($arr[2] ne $chr){
		return "bad";
	}
	my $len=length($arr[9]);
	if($arr[3]+$len < $start){
		return "bad";
	}
	if($arr[3] > $end){
		return "bad";
	}
	return "good";
}
