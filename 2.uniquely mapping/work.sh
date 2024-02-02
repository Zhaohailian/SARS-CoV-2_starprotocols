perl step1.prepare_reads_fragment.pl rm_FalsePosi.sam
hisat2 -p 12 -k 2 -x human_virus.pan.fa -U rm_FalsePosi.fq -S rm_FalsePosi.reMapHisat2.sam --un rm_FalsePosi.unMapHisat2.fq 2> rm_FalsePosi.hisat2.log
bowtie2 -un rm_FalsePosi.unMapBowtie.fq --sensitive -N 1 -L 17 --end-to-end -p 12 -x human_virus.pan.fa -U rm_FalsePosi.fq -S rm_FalsePosi.reMapBowtie.sam 2> rm_FalsePosi.reMapBowtie2.log
STAR --runMode alignReads --genomeDir Pangenome --readFilesIn rm_FalsePosi.fq --outFileNamePrefix rm_FalsePosi.toGenome_ --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.9 --outFilterScoreMinOverLread 0.9 --outSAMattributes All --runThreadN 16
perl step2.separatBothUniq_and_Other.requireLocitoo.pl rm_FalsePosi.fq rm_FalsePosi.sam rm_FalsePosi.reMapHisat2.sam rm_FalsePosi.reMapBowtie.sam rm_FalsePosi.toGenome_Aligned.out.sam
