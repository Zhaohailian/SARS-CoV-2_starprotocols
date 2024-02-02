perl find_targets_related_to_a_region_plus.pl intergene.rm_FalsePosi.bothUniq.rm_Random.sam  NC_045512.2 0 30000 > virus.TargetAndSource.sam
perl extract_junction.pl virus.TargetAndSource.sam
python junction_fragment_region.py out1.arms.info 50 > fragment50nt.bed
sort -k1,1 -k2,2n fragment50nt.bed |awk -F '\t' 'BEGIN{OFS="\t"{if($2>=0) print $0 }' > fragment50nt.sort.bed 
grep -E 'PC5UTR|PCCDS|PC3UTR' chlSab2.gene_element.bed |bedtools intersect -a fragment50nt.sort.bed -b - -s -u > fragment50nt.exon.bed
Piranha -s -o peaks.bed -p 0.05 -b 20 -a 0.96 fragment50nt.exon.bed
