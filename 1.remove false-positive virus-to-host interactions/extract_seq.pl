#!/usr/bin/perl
use Array::Compare;

die "perl $0 interaction.sam" if(@ARGV != 1);
my $sam=shift;

open(OUTA,">part_A.fa") || die;
open(OUTB,">part_B.fa") || die;

my $count;
open(SM,$sam) || die;
while(my $frag_a=<SM>){
        if($frag_a=~/^@/){
                next;
        }
        else{
                my $frag_b=<SM>;
		my @sub_a=split/\s+/,$frag_a;
		my @sub_b=split/\s+/,$frag_b;
                my @id_a=split/_/,$sub_a[0];
                my @id_b=split/_/,$sub_b[0];

                if($id_a[0]."_".$id_a[1] ne $id_b[0]."_".$id_b[1]){
                        die "wrong format\n";
                }
                else{#same read name
                        my $chr_a=$sub_a[2];
                        my $loci_a=$sub_a[3];
                        my $strand_a=$sub_a[1];
                        my $cigar_a=$sub_a[5];
                        my $seq_a=$sub_a[9];

                        my $chr_b=$sub_b[2];
                        my $loci_b=$sub_b[3];
                        my $strand_b=$sub_b[1];
			my $cigar_b=$sub_b[5];
                        my $seq_b=$sub_b[9];

			if($sub_a[2] !~ /NC_045512.2/ and $sub_b[2] !~ /NC_045512.2/){	#at least one arm in NC
				next;
			}

			my ($start_match_a,$match_len_a)=get_start_and_len($cigar_a);
			my ($start_match_b,$match_len_b)=get_start_and_len($cigar_b);

			if($start_match_a eq "NA" or $start_match_b eq "NA"){
				#warn $frag_a,$frag_b;
				next;
			}

			print OUTA ">",$id_a[0]."_".$id_a[1],"\t",$cigar_a,"\n";
			print OUTA substr($seq_a,$start_match_a,$match_len_a),"\n";
			print OUTB ">",$id_b[0]."_".$id_b[1],"\t",$cigar_b,"\n";
			print OUTB substr($seq_b,$start_match_b,$match_len_b),"\n";
		}
	}
}

sub get_start_and_len{
        my $raw_cigar=shift;
        my $start;
        my $only_block_len=0;
        my $only_block_start=0;
	my $cigar=$raw_cigar;
        while($cigar=~/(\d+)(\w)/g){
		my $num=$1;
		my $class=$2;
		if($class eq "S" or $class eq "H"){
			$only_block_start+=$1;
		}
		else{
			last;
		}
	}
	my $cigar=$raw_cigar;
	while($cigar=~/(\d+)(\w)/g){
                my $num=$1;
                my $class=$2;
                if($class eq "S" or $class eq "H"){
                #        $only_block_start+=$1;
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
        if($only_block_len < 10){
                return ("NA","NA");
        }
        else{
                return ($only_block_start,$only_block_len)
        }
}
	








