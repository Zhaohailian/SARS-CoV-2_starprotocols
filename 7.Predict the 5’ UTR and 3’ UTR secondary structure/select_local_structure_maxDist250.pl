#!/usr/bin/perl
die "perl $0 3UTR.ct  update.mergeBwa.virus.resolution2.vcrt.format.matrix\n" if(@ARGV != 2);
my $ct_file=shift;
my $RICseq_count_matrix_for_correlation=shift;

#my $extended_loop_start=29545;   ### 3UTR
#my $extended_loop_end=29870;
#my $loop_start=29545;
#my $loop_end=29870;

my $extended_loop_start=0;  ### 5UTR
my $extended_loop_end=400;
my $loop_start=0;
my $loop_end=400;

my $ruler_used=2; 
my $print_paired_pixels=1;
my $print_pvalue_eachStructure=1;
my $output_sequence_id="NC_045512.2";

#read RIC-seq count matrix
my %RICseq_matrix;
open(MT,$RICseq_count_matrix_for_correlation) || die;
my $head=<MT>;
chomp $head;
my @arr_head=split/\s+/,$head;
while(my $line=<MT>){
        chomp $line;
        my @sub=split/\s+/,$line;

        foreach (1..$#sub){
                if($arr_head[$_] <= $sub[0]){
                        next;
                }
                $RICseq_matrix{$sub[0]}{$arr_head[$_]}=$sub[$_];
        }
}

#finished RIC-seq count matrix

my ($index,$energy,$pvalue,$ctinfo)=select_by_correlation($extended_loop_start,$extended_loop_end,$ct_file,$loop_start,$loop_end);
my @ct_lines=split/\n/,$ctinfo;
shift @ct_lines;
foreach my $l (@ct_lines){
     	$l=~s/^\s+//;
      	my @sub=split/\s+/,$l;
     	if($#sub != 5){
        	next;
     	}
  	if($sub[4]){
      		my $real_loci_left=$extended_loop_start+$sub[0];
    		my $real_loci_right=$extended_loop_start+$sub[4];
   		my @pair=($real_loci_left,$real_loci_right);
     		@pair=sort {$a<=>$b} @pair;
 		$dynamic_constrains{$pair[0]."\t".$pair[1]}=$pair[0];
  	}
}

sub select_by_correlation{
	my $start=shift;
	my $end=shift;
	my $ct_file=shift;
	my $this_loop_start=shift;
	my $this_loop_end=shift;
	my $slop_due_to_RICseq_resolution=0;
	my $least_distance_to_diagonal=5;
	my %structure;
	my $id_of_structure;
	$/="ENERGY";
	open(CT,$ct_file) || die;
	<CT>;
	while(my $block=<CT>){
		chomp $block;
		$id_of_structure++;
		my %paired;
		my @lines=split/\n/,$block;
		my $first=shift @lines;
		$first=~s/^\s+//;
		my $energy=(split/\s+/,$first)[1];
		foreach my $l (@lines){
			$l=~s/^\s+//;
			my @sub=split/\s+/,$l;
			if($sub[4]){#paired
				$paired{$sub[0]}{$sub[4]}=1;
			}
			else{#single strand
			}
		}
		
		#do t-test
		my %paired_RICseq_signal;
		foreach my $i ($start..$end){
			my $win_i=int($i/$ruler_used)*$ruler_used;
			foreach my $j ($start..$end){
				my $win_j=int($j/$ruler_used)*$ruler_used;
				if($win_j <= $win_i){
					next;
				}
				else{
					my $relative_loci_i=$i-$start+1;
					my $relative_loci_j=$j-$start+1;
					if($paired{$relative_loci_i}{$relative_loci_j}){
						$paired_RICseq_signal{$win_i."\t".$win_j}=1;
					}
				}
			}
		}
		my @paired;
		my @other;
		my %unique_windows;
		foreach my $i ($start..$end-1){
			my $win_i=int($i/$ruler_used)*$ruler_used;
			foreach my $j ($start..$end-1){
				my $win_j=int($j/$ruler_used)*$ruler_used;
				if($win_j <= $win_i+$least_distance_to_diagonal){
					next;
				}
				if($unique_windows{$win_i."\t".$win_j}){	#already assigned
					next;
				}
				if($resolved_bases{$i} and $resolved_bases{$j}){	#keep focus on newly added region
					next;
				}
				if($i < $this_loop_start or $i > $this_loop_end or $j < $this_loop_start or $j > $this_loop_end){	#keep focus on newly added by this loop
					next;
				}
				$unique_windows{$win_i."\t".$win_j}=1;
				my $signal=$RICseq_matrix{$win_i}{$win_j};
				if(defined($signal)){
				}
				else{
					warn $ct_file,"\t",$i,"\t",$j,"\t",$win_i,"\t",$win_j,"\tbug\n";
					die;
				}
				if(exists $paired_RICseq_signal{$win_i."\t".$win_j}){
					push (@paired,$signal);
				}
				else{
					my $distance_to_paired=2*$slop_due_to_RICseq_resolution+1;
					foreach my $paired_pixels (keys %paired_RICseq_signal){
						my ($paired_win_i,$paired_win_j)=split/\s+/,$paired_pixels;
						my $tmp_dis;
						if(abs($win_i-$paired_win_i) > abs($win_j-$paired_win_j)){
							$tmp_dis=abs($win_i-$paired_win_i);
						}
						else{
							$tmp_dis=abs($win_j-$paired_win_j);
						}
						if($tmp_dis < $distance_to_paired){
							$distance_to_paired=$tmp_dis;
						}
					}
					if($distance_to_paired <= $slop_due_to_RICseq_resolution){
						push (@paired,$signal);
					}
					else{
						push (@other,$signal);
					}
				}
			}
		}
		
		open(TMPP,">tmp.paired.RICseq.signal.list") || die;
		open(TMPO,">tmp.other.RICseq.signal.list") || die;
		foreach (@paired){
			print TMPP $_,"\n";
		}
		foreach (@other){
			print TMPO $_,"\n";
		}
		close TMPP;
		close TMPO;

		if(@paired < 3 or @other < 3){
			next;
		}

		my $test_result=`Rscript 0.Rscript_ttest.r`;
		$test_result=~s/\s+$//;
		my ($pvalue,$mean_paired,$median_paired,$mean_other,$median_other)=split/\s+/,$test_result;

		my $paired_sum=$mean_paired*($#paired+1);
		my $other_sum=$mean_other*($#other+1);

		$structure{$id_of_structure}{"energy"}=$energy;
		$structure{$id_of_structure}{"sum"}=$paired_sum;
		$structure{$id_of_structure}{"pvalue"}=$pvalue;
		$structure{$id_of_structure}{"ct"}=$block;
		if($print_pvalue_eachStructure){
			print ">",$id_of_structure,"\t",$#paired+1,"\t",$paired_sum,"\t",$mean_paired,"\t",$median_paired,"\t",$#other+1,"\t",$other_sum,"\t",$mean_other,"\t",$median_other,"\t",$pvalue,"\ttestPvalue\n";
		}
	}
	$/="\n";

	if(%structure){
		#my @sorted_strcuture=sort {$structure{$a}{"pvalue"} <=> $structure{$b}{"pvalue"}} keys %structure;	#min pvalue
		my @sorted_strcuture=sort {$structure{$b}{"sum"} <=> $structure{$a}{"sum"}} keys %structure;		#max sum
		#return ($sorted_strcuture[0],$structure{$sorted_strcuture[0]}{"energy"},$structure{$sorted_strcuture[0]}{"pvalue"},$structure{$sorted_strcuture[0]}{"ct"});

		my $optimal_value=$structure{$sorted_strcuture[0]}{"sum"};

		my @candidate_structure;
		foreach (@sorted_strcuture){
			if($structure{$_}{"sum"} >= $optimal_value){
				push (@candidate_structure,$_);
			}
		}
		
		my @sorted_again_structure=sort {$structure{$a}{"energy"} <=> $structure{$b}{"energy"}} @candidate_structure;	#lowest energy
		#my @sorted_again_structure=(1);	#test
		return ($sorted_again_structure[0],$structure{$sorted_again_structure[0]}{"energy"},$structure{$sorted_again_structure[0]}{"pvalue"},$structure{$sorted_again_structure[0]}{"ct"});
	}
	else{
		return ("NA","NA","NA","NA");
	}

}

	

