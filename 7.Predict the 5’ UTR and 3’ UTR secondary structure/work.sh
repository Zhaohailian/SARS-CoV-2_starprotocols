perl sam_to_loops.pl virus_chimeric_reads.sam | sort -k3,3 -k7,7 > virus.sort.list
java -jar juicebox_tools.jar pre -r 1,2,5,10,25,50,100 virus.sort.list virus.hic hs.size
java -jar juicebox_tools.jar dump observed VC_SQRT virus.hic hs1 hs1 BP 2 virus.2nt.vcrt.matrix
perl format_hic_to_matrix.pl virus.2nt.vcrt.matrix > virus.2nt.vcrt.format.matrix
Fold 3UTR.fa 3UTR.ct --maxdistance 250
perl select_local_structure.pl 3UTR.ct virus.2nt.vcrt.format.matrix > log.3UTR.txt
python3 select_structure.py log.3UTR.txt 3UTR.ct > best_structure.ct
