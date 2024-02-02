perl from_sam_to_pair_reads_bed.pairs.pl intergene.rm_FalsePosi.bothUniq.sam
bedtools pairtopair -a intergene.rm_FalsePosi.bothUniq.bedpair -b random.intergene.rm_FalsePosi.bothUniq.bedpair -type notboth -f 0.05 > intergene.rm_FalsePosi.bothUniq.rm_Random.bedpair
python obtain_sam_file.py intergene.rm_FalsePosi.bothUniq.rm_Random.bed intergene.rm_FalsePosi.bothUniq.sam > intergene.rm_FalsePosi.bothUniq.rm_Random.sam

