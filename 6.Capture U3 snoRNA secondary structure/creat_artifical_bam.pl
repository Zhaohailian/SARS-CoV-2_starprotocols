#!/usr/bin/perl
die "perl $0 virus.onlyPart.add_readLoci.withID.list out4.onlyPart.loops_of_enriched_pixels.resolution2.score0.01.bedpe ref.fa\n" if(@ARGV != 3);
my $read_arm_info_table=shift;
my $cluster_loops_table=shift;
my $ref_fasta=shift;

my $refseq;
open(REF,$ref_fasta) || die;
while(my $line=<REF>){
	chomp $line;
	if($line=~/>/){
		next;
	}
	$refseq.=$line;
}

my %read_arms;
open(RAIT,$read_arm_info_table) || die;
while(my $line=<RAIT>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$read_arms{$sub[0]}=$line;
}

open(OUTP,">U3.read_of_cluster.sam") || die;
open(CLT,$cluster_loops_table) || die;
<CLT>;
while(my $line=<CLT>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($sub[6] !~ /ManualLoops/){
		next;
	}
	my $id=$sub[6];
	foreach (keys %read_arms){
		my @info=split/\s+/,$read_arms{$_};
		if($info[2] > $sub[1] and $info[2] < $sub[2] and $info[5] > $sub[4] and $info[5] < $sub[5]){
			if($info[-1] =~ /Chimeric/){
				next;
				if(($info[11]-$info[9]-1)>=0){
					print OUTC $info[-1],"\t0\tNR_006880.1\t",$info[8],"\t255\t",$info[9]-$info[8]+1,"M";
					print OUTC $info[11]-$info[9]-1,"N",$info[12]-$info[11]+1,"M\t";
					print OUTC "*\t0\t0\t";
					print OUTC substr($refseq,$info[8]-1,$info[9]-$info[8]+1);
					print OUTC substr($refseq,$info[11]-1,$info[12]-$info[11]+1);
					my $len=$info[9]-$info[8]+1+$info[12]-$info[11]+1;
					my $qual="F" x $len;
					print OUTC "\t$qual\t";
					print OUTC "RG:Z:$id\n";
				}
			}
			elsif($info[-1] =~ /Part/){
				if(($info[11]-$info[9]-1)>=0){
	                                print OUTP $info[-1],"\t0\tNR_006880.1\t",$info[8],"\t255\t",$info[9]-$info[8]+1,"M";
					print OUTP $info[11]-$info[9]-1,"N",$info[12]-$info[11]+1,"M\t";
					print OUTP "*\t0\t0\t";
					print OUTP substr($refseq,$info[8]-1,$info[9]-$info[8]+1);
					print OUTP substr($refseq,$info[11]-1,$info[12]-$info[11]+1);
					my $len=$info[9]-$info[8]+1+$info[12]-$info[11]+1;
					my $qual="F" x $len;
					print OUTP "\t$qual\t";
					print OUTP "RG:Z:$id\n";
				}
			}
			else{
				die;
			}
		}
	}
	close OUTA;
}






