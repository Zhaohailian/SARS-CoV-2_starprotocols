perl extract_seq.pl intergene.sam
blastn -query part_A.fa -task megablast -db virus.fa -out part_A.out -evalue 1 -word_size 5 -outfmt 7 -num_threads 16
blastn -query part_B.fa -task megablast -db virus.fa -out part_B.out -evalue 1 -word_size 5 -outfmt 7 -num_threads 16
perl filter_according_to_BLASTN.pl part_A.out part_B.out intergene.sam 0.8 > intergene.rm_FalsePosi.sam

