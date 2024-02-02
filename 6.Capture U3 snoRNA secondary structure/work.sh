perl sam_to_bedpe.pl U3_chimeric_reads.sam 0.01 U3_chimeric_reads.list U3_chimeric_reads.2nt.bedpe
perl cluster_pixels.pl U3_chimeric_reads.2nt.bedpe > U3_chimeric_reads.cluster.bedpe
perl creat_artifical_bam.pl U3_chimeric_reads.list U3_chimeric_reads.cluster.bedpe U3.fasta
awk -F '\t' 'BEIGN{OFS="\t"}{if($12=="RG:Z:ManualLoops_1") print $0 }' U3.read_of_cluster.sam > U3.read_of_cluster1.sam
samtools view -t chr.sizes -bSo U3.read_of_cluster1.bam U3.read_of_cluster1.sam
samtools sort -o U3.read_of_cluster1.sort.bam U3.read_of_cluster1.bam
samtools index U3.read_of_cluster1.sort.bam
