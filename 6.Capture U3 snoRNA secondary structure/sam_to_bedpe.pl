#!/usr/bin/perl
die "perl $0 pcp_rep1_interaction.sam want_chr want_start want_end> xxx.bedpe\n" if(@ARGV != 4);
my $contact_sam=shift;
my $want_chr="NR_006880.1";
my $resolution=2;   ### bedpe resolution
my $length=217;   ### length of U3 snoRNA
my $score_cutoff=shift;
my $minimal_loop_len=4*$resolution;  ### minimal loop length
my $output_arm_list=shift;  ## output1: gaps and arms of part type chimeric reads
my $output_pixel=shift; ## output2: pixel bedpe with connection score than score_cutoff

my $sampleName=$contact_sam;
$sampleName=~s/^.+\///;
$sampleName=~s/_interaction.Chimeric.HQ.sam//;

my $start_win_index=0;
my $end_win_index=int($length/$resolution);
 
my $num_of_unique;
my %unique;
my %pairwise_unique;
my %pairwise;
my %coverage;

my $read_id;
my $pixel_id;
my %unique;
my %gap_len;
open(OUTPUT1, ">$output_arm_list") or die "Cannot open output file\n!";
open(OUTPUT2, ">$output_pixel") or die "Cannot open output file\n!";

open(CS,$contact_sam) || die;
while(my $frag_a=<CS>){
        if($frag_a=~/^@/){
                next;
        }
        my $frag_b=<CS>;
	my @sub_a=split/\s+/,$frag_a;
	my @sub_b=split/\s+/,$frag_b;
	my @id_info_a=split/_/,$sub_a[0];
	my @id_info_b=split/_/,$sub_b[0];

     	if($id_info_a[0]."_".$id_info_a[1] ne $id_info_b[0]."_".$id_info_b[1]){
        	die "wrong format\n";
      	}
     	else{#same read name
        	my $chr_a=$sub_a[2];
		my ($loci_a,$end_a)=get_map_start_end($frag_a);
		$sub_a[1] = $sub_a[1] > 255 ? $sub_a[1]-256 : $sub_a[1];

  		my $chr_b=$sub_b[2];
		my ($loci_b,$end_b)=get_map_start_end($frag_b);
		$sub_b[1] = $sub_b[1] > 255 ? $sub_b[1]-256 : $sub_b[1];

		if($id_info_a[2] ne $id_info_b[2] or $id_info_a[2] ne "Plus"){	#same reads; same RNA; same strand
			next;
		}		
		if($chr_a ne $want_chr or $chr_b ne $want_chr){
			next;
		}
                if($sub_a[0] =~ /^AlignPair/){
			next;
                }
                elsif($sub_a[0] =~ /^Part/){
                        if($end_a >= $loci_b){	#impossible situation
                                print $frag_a,$frag_b;
                                die "this Part-pair is incorrect\n";
                        }
			my $tmp_gap_len=$loci_b-$end_a-1;
			$gap_len{$tmp_gap_len}++;
			$read_id++;
			my $line1;
			$line1=join("\t", $read_id,$want_chr,$end_a,$end_a+1,$want_chr,$loci_b-1,$loci_b,$want_chr,$loci_a,$end_a,$want_chr,$loci_b,$end_b,$sampleName,"Part\tMapped",$id_info_a[0]."_".$id_info_a[1]);
			print OUTPUT1 $line1,"\n";
			
			$arms_loci{$read_id}=$line1;
			my $win_a=int($end_a/$resolution)*$resolution;
			my $win_b=int(($loci_b-1)/$resolution)*$resolution;
			if($win_a eq $win_b){
				$coverage{$win_a}++;
				$pairwise{$win_a}{$win_b}++;
			}
			else{
				$coverage{$win_a}++;
				$coverage{$win_b}++;
				$pairwise{$win_a}{$win_b}++;
				$pairwise{$win_b}{$win_a}++;
			}
			my @arm_info=split/\s+/,$arms_loci{$read_id};
			my $arms=join"\t",@arm_info[7..15];
			if($unique{$arms}){
				next;
			}
			else{
				$unique{$arms}=1;
				$num_of_unique++;
				if($win_a eq $win_b){
					$pairwise_unique{$win_a}{$win_b}++;
				}
				else{
					$pairwise_unique{$win_a}{$win_b}++;
					$pairwise_unique{$win_b}{$win_a}++;
				}
			}
		}
                elsif($sub_a[0] =~ /^Chimeric/){
			next;
		}
		else{
			print $frag_a,$frag_b;
			die "wrong reads class\n";
		}
	}
}

warn $num_of_unique," unique pair tags\n";

my $header="#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2";
print OUTPUT2 $header,"\n";

foreach my $i ($start_win_index..$end_win_index){
	my $win_a=$i*$resolution;
	foreach my $j ($start_win_index..$end_win_index){
		if($j < $i){
			next;
		}
		my $win_b=$j*$resolution;
		if(!$pairwise{$win_a}{$win_b}){
			next;
		}
		my $depth_a=$coverage{$win_a};
		my $depth_b=$coverage{$win_b};
		my $sqrt_depth=sqrt($depth_a*$depth_b);
		my $score=$pairwise{$win_a}{$win_b}/$sqrt_depth;
		if($win_b - $win_a < 0){
			die;
		}
		if($win_b - $win_a < $minimal_loop_len){	#minmal loop len
			next;
		}
		if($pairwise_unique{$win_a}{$win_b} < 2){	#unique reads
			next;
		}
		if($score < $score_cutoff){	#score
			next;
		}
		$pixel_id++;
		#print "hs1\t$win_a\t",$win_a+$resolution,"\t";
		#print "hs1\t$win_b\t",$win_b+$resolution,"\t";
		#print "EnrichPixel_$pixel_id\t255\t+\t+\n";
		my $out_line="hs1\t$win_a\t".($win_a+$resolution)."\ths1\t".$win_b."\t".($win_b+$resolution)."\tEnrichPixel_".$pixel_id."\t255\t+\t+";
		print OUTPUT2 $out_line,"\n";

	}
}


sub get_map_start_end{
        my $frag=shift;
        my @sub=split/\s+/,$frag;
        my $cigar=$sub[5];
        my $chr_start=$sub[3];
        my $chr_end=$sub[3];
        while($cigar=~/(\d+)(\w)/g){
                my $tmp_len=$1;
                my $tmp_content=$2;
                if($tmp_content eq "M"){
                        $chr_end+=$tmp_len;
                }
                elsif($tmp_content eq "I"){
                }
                elsif($tmp_content eq "D"){
                        $chr_end+=$tmp_len;
                }
                elsif($tmp_content eq "S" or $tmp_content eq "H"){
                }
                elsif($tmp_content eq "N"){
                        warn "No gap exists in chimeric segments\n";
                        warn $frag;
                        die;
                }
        }
        return ($chr_start,$chr_end-1);
}

