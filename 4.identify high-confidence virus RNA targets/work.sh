perl MonteCarlo_simulation_use_fragment.pl whole_gene_region.bed 100000 0.05 20 intergene_simulation intergene.rm_FalsePosi.bothUniq.rm_Random.sam > run.sh
bash run.sh
awk -F '\t' 'BEGIN{OFS="\t"}{if(($1=="NC_045512.2")&&($13>=2)) print $7,$8,$9,$10,$11,$12;else if(($7=="NC_045512.2")&&($13>=2)) print $1,$2,$3,$4,$5,$6 }' intergene_simulation.significant.interMolecular.interaction.list > virus_target.txt 
