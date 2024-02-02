#!/usr/bin/perl
die "perl $0 out1.virus_to_exons.sam length_of_interaction_region\n" if(@ARGV != 1);
my $chimeric_reads_sam=shift;
#my $target_region_len=shift;

#open(VTAG,">out2.tag_on_virs.bed") || die;
#open(GTAG,">out2.tag_on_genome.bed") || die;
open(INFO,">out1.arms.info") || die;	#help to assemble chimeric fragments

open(SM,$chimeric_reads_sam) || die;
while(my $frag_a=<SM>){
        if($frag_a=~/^@/){
                next;
        }
        else{
		my $frag_b=<SM>;
                my @sub_a=split/\s+/,$frag_a;
                my @sub_b=split/\s+/,$frag_b;
                my @id_a_info=split/_/,$sub_a[0];
                my @id_b_info=split/_/,$sub_b[0];
                my $sample_a=shift @id_a_info;
                my $sample_b=shift @id_b_info;
                my $strand_a=$id_a_info[2];
                my $strand_b=$id_a_info[2];

                my $chr_a=$sub_a[2];
                my $loci_a=$sub_a[3];
                my $cigar_a=$sub_a[5];
                $cigar_a=~/(\d+)M/;
                my $match_a=$1;

                my $chr_b=$sub_b[2];
                my $loci_b=$sub_b[3];
                my $cigar_b=$sub_b[5];
                $cigar_b=~/(\d+)M/;
                my $match_b=$1;


		my $readId=$sample_a."_".$id_a_info[0]."_".$id_a_info[1];
	
		if($id_a_info[1] ne $id_b_info[1]){
			die;
		}
		if($id_a_info[0] !~ /Chimeric/ or $id_b_info[0] !~ /Chimeric/){
			next;
		}
		if($sub_a[2] =~ /NC_045512.2/ and $sub_b[2] !~ /NC_045512.2/){
			if($id_a_info[2] ne "Plus"){	#stand specific
				next;
			}
			#seq on virus
			if($sub_a[1] == 256 or $sub_a[1] == 0){
				my $junction=$loci_a+$match_a-1;	#get seq at the left of the junction
				#print VTAG "NC_045512.2\t",$junction-$target_region_len,"\t",$junction,"\t",$readId,"\t255\t+\n";
				print INFO "NC_045512.2\t",$junction,"\t";
				print INFO "+\t";
			}
			elsif($sub_a[1] == 272 or $sub_a[1] == 16){
				my $junction=$loci_a;		#get seq at the right of the junction
				#print VTAG "NC_045512.2\t",$junction-1,"\t",$junction-1+$target_region_len,"\t",$readId,"\t255\t+\n";
				print INFO "NC_045512.2\t",$junction,"\t";
				print INFO "-\t";
			}
			else{
				die;
			}

			#seq on genome
			my $RNA_strand=$id_b_info[2];
			if($sub_b[1] == 256 or $sub_b[1] == 0){
				if($RNA_strand eq "Plus"){
					my $junction=$loci_b;	#get seq at the right of the junction; no need to complementary & reverse
					#print GTAG $sub_b[2],"\t",$junction-1,"\t",$junction-1+$target_region_len,"\t",$readId,"\t255\t+\n";
					print INFO $sub_b[2],"\t",$junction,"\t";
					print INFO "+\t$RNA_strand\t";
				}
				elsif($RNA_strand eq "Minus"){
					my $junction=$loci_b;	#get seq at the right of the junction; need to complementary & reverse
					#print GTAG $sub_b[2],"\t",$junction-1,"\t",$junction-1+$target_region_len,"\t",$readId,"\t255\t-\n";
					print INFO $sub_b[2],"\t",$junction,"\t";
					print INFO "+\t$RNA_strand\t";
				}
				else{
					die;
				}
			}
			elsif($sub_b[1] == 272 or $sub_b[1] == 16){
				if($RNA_strand eq "Plus"){
					my $junction=$loci_b+$match_b-1; #get seq at the left of the junction; no need to complementary & reverse
					#print GTAG $sub_b[2],"\t",$junction-$target_region_len,"\t",$junction,"\t",$readId,"\t255\t+\n";
					print INFO $sub_b[2],"\t",$junction,"\t";
					print INFO "-\t$RNA_strand\t";
				}
				elsif($RNA_strand eq "Minus"){
					my $junction=$loci_b+$match_b-1; #get seq at the left of the junction; need to complementary & reverse
					#print GTAG $sub_b[2],"\t",$junction-$target_region_len,"\t",$junction,"\t",$readId,"\t255\t-\n";
					print INFO $sub_b[2],"\t",$junction,"\t";
					print INFO "-\t$RNA_strand\t";
				}
				else{
					die;
				}
			}
			else{
				die;
			}
			print INFO "VirusHeadGenomeTail\t$readId\n";
		}
		elsif($sub_a[2] !~ /NC_045512.2/ and $sub_b[2] =~ /NC_045512.2/){
			if($id_b_info[2] ne "Plus"){	#strand specific
				next;	
			}
			#seq on virus
			if($sub_b[1] == 256 or $sub_b[1] == 0){
				my $junction=$loci_b;		#get seq at the right of the junction
				#print VTAG "NC_045512.2\t",$junction-1,"\t",$junction-1+$target_region_len,"\t",$readId,"\t255\t+\n";
				print INFO "NC_045512.2\t",$junction,"\t";
				print INFO "+\t";
			}
			elsif($sub_b[1] == 272 or $sub_b[1] == 16){
				my $junction=$loci_b+$match_b-1;	#get seq at the left of the junction
				#print VTAG "NC_045512.2\t",$junction-$target_region_len,"\t",$junction,"\t",$readId,"\t255\t+\n";
				print INFO "NC_045512.2\t",$junction,"\t";
				print INFO "-\t";
			}
			else{
				die;
			}
	
			#seq on genome
			my $RNA_strand=$id_a_info[2];
			if($sub_a[1] == 256 or $sub_a[1] == 0){
				if($RNA_strand eq "Plus"){
					my $junction=$loci_a+$match_a-1; #get seq at the left of the junction; no need to complementary & reverse
					#print GTAG $sub_a[2],"\t",$junction-$target_region_len,"\t",$junction,"\t",$readId,"\t255\t+\n";
					print INFO $sub_a[2],"\t",$junction,"\t";
					print INFO "+\t$RNA_strand\t";
				}
				elsif($RNA_strand eq "Minus"){
					my $junction=$loci_a+$match_a-1; #get seq at the left of the junction; need to complementary & reverse
					#print GTAG $sub_a[2],"\t",$junction-$target_region_len,"\t",$junction,"\t",$readId,"\t255\t-\n";
					print INFO $sub_a[2],"\t",$junction,"\t";
					print INFO "+\t$RNA_strand\t";
				}
				else{
					die;
				}
			}
			elsif($sub_a[1] == 272 or $sub_a[1] == 16){
				if($RNA_strand eq "Plus"){
					my $junction=$loci_a;	#get seq at the right of the junction; no need to complementary & reverse
					#print GTAG $sub_a[2],"\t",$junction-1,"\t",$junction-1+$target_region_len,"\t",$readId,"\t255\t+\n";
					print INFO $sub_a[2],"\t",$junction,"\t";
					print INFO "-\t$RNA_strand\t";
				}
				elsif($RNA_strand eq "Minus"){
					my $junction=$loci_a;	#get seq at the right of the junction; need to complementary & reverse
					#print GTAG $sub_a[2],"\t",$junction-1,"\t",$junction-1+$target_region_len,"\t",$readId,"\t255\t-\n";
					print INFO $sub_a[2],"\t",$junction,"\t";
					print INFO "-\t$RNA_strand\t";
				}
				else{
					die;
				}
			}
			else{
				die;
			}
			print INFO "GenomeHeadVirusTail\t$readId\n";
		}
		else{
			die;
		}
	}
}





	
