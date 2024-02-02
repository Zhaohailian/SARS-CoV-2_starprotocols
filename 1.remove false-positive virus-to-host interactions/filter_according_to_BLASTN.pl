#!/usr/bin/perl
die "perl $0 part_A_match_seq_to_rRNA.out part_B_match_seq_to_rRNA.out GM12878_rep1_interaction.without_rRNA.unique.rmPoly.sam > GM12878_rep1_interaction.without_rRNA.unique.rmPoly.rmRNA.sam ratio_cutoff\n" if(@ARGV != 4);
my $part_a_rRNA=shift;
my $part_b_rRNA=shift;
my $in_sam=shift;
my $ratio_cutoff=shift;

my %bad_a;
open(RNAA,$part_a_rRNA) || die;
while(my $line=<RNAA>){
	if($line=~/^#/){
		next;
	}
	else{
		my @sub=split/\s+/,$line;
		if(exists $bad_a{$sub[0]}){
			next;
		}
		else{
			$bad_a{$sub[0]}=$sub[3];
		}
	}
}

my %bad_b;
open(RNAB,$part_b_rRNA) || die;
while(my $line=<RNAB>){
        if($line=~/^#/){
                next;
        }
        else{
                my @sub=split/\s+/,$line;
		if(exists $bad_b{$sub[0]}){
			next;
		}
		else{
			$bad_b{$sub[0]}=$sub[3];
		}
        }
}


open(BAD,">FalsePositive_linked.sam") || die;
open(SM,$in_sam) || die;
while(my $frag_a=<SM>){
        if($frag_a=~/^@/){
		print $frag_a;
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
			my @sub_a=split/\s+/,$frag_a;
                        my $cigar_a=$sub_a[5];
                        my $match_a=cigar_to_match_len($cigar_a);
			my $ratio_a;
			if(exists $bad_a{$id_a[0]."_".$id_a[1]}){
				$ratio_a=$bad_a{$id_a[0]."_".$id_a[1]}/$match_a;
			}
			else{
				$ratio_a=0;
			}

                        my @sub_b=split/\s+/,$frag_b;
                        my $cigar_b=$sub_b[5];
                        my $match_b=cigar_to_match_len($cigar_b);
			my $ratio_b;
			if(exists $bad_b{$id_b[0]."_".$id_b[1]}){
				$ratio_b=$bad_b{$id_b[0]."_".$id_b[1]}/$match_b;
			}
			else{
				$ratio_b=0;
			}
			
			if($ratio_a >= $ratio_cutoff and $ratio_b >= $ratio_cutoff and ($sub_a[2] =~ /NC_045512.2/ or $sub_b[2] =~ /NC_045512.2/)){
				print BAD $frag_a,$frag_b;
				next;
			}
			else{
				print $frag_a,$frag_b;
			}
		}
	}
}
	
	
sub cigar_to_match_len{
        my $cigar=shift;
        my $only_block_len;
        my $only_block_start;
        while($cigar=~/(\d+)(\w)/g){
                my $num=$1;
                my $class=$2;
                if($class eq "S" or $class eq "H"){
                        $only_block_start+=$1;
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
	return $only_block_len;
}

