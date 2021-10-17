#!/bin/bash
#SBATCH --job-name=variant_calling_merged_gatk2
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=4000M
#SBATCH -p long
#SBATCH --output=variant_calling_merged_gatk2.out
#SBATCH --error=variant_calling_merged_gatk2.error
#SBATCH --mail-user=m.corteslopez@imb-mainz.de
#SBATCH --mail-type=ALL

while read line; 
	do
		BC=$line
		gatk --java-options "-Xmx4g" HaplotypeCaller -R ref/minigene.fa -I merged_aln/${BC}.bam -O merged_vcfs/${BC}.vcf.gz
done< total_merged.bc
