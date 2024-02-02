#!/usr/bin/perl
die "perl virus.2nt.vcrt.matrix > virus.2nt.vcrt.format.matrix \n" if(@ARGV != 1);
my $in_matrix=shift;
my $chr_size=14951;       
my $ruler=2;


my %matrix;
open(IN,$in_matrix) || die;
while(my $line=<IN>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$matrix{$sub[0]}{$sub[1]}=$sub[2];
	$matrix{$sub[1]}{$sub[0]}=$sub[2];
}

my @all;
foreach (0..$chr_size){
	push (@all,$ruler*$_);
}
	

print "loci\t";

foreach (@all){
	print $_,"\t";
}
print "\n";

foreach my $row (@all){
	print $row,"\t";
	foreach my $col (@all){
		print $matrix{$row}{$col}+0,"\t";
	}
	print "\n";
}
