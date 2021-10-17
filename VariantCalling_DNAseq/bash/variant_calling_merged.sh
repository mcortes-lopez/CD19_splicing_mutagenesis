#!/bin/bash
#SBATCH --job-name=variant_calling_merged
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=4000M
#SBATCH -p long
#SBATCH --output=variant_calling_merged.out
#SBATCH --error=variant_calling_merged.error
#SBATCH --mail-user=m.corteslopez@imb-mainz.de
#SBATCH --mail-type=ALL


source activate pacbio
module load samtools/1.3.1
module load jdk/1.8.0_181


totallines=$(wc -l total_merged.bc | awk '{print $1}')

i=0
while [ $i -le $totallines ]
	do
	awk -v st="$i" 'NR>=st && NR< st+100{print $0}' total_merged.bc > tmp100lines
	while read line; 
		do
			BC=$line
			gunzip merged_fastq/CD19_ccs_$BC.fastq.gz
			blasr merged_fastq/CD19_ccs_$BC.fastq ref/minigene.fa --bam --out $BC.bam
			gzip merged_fastq/CD19_ccs_$BC.fastq
			samtools sort $BC.bam > $BC.sorted.bam
			samtools view -H $BC.sorted.bam | awk -F '\t' '{if($0~/^@RG/){print $0"\tSM:sample1"}else{print $0}}' | samtools reheader - $BC.sorted.bam > merged_aln/$BC.bam
			samtools index merged_aln/$BC.bam
			rm $BC.sorted.bam $BC.bam &

#			gatk --java-options "-Xmx4g" HaplotypeCaller -R ref/minigene.fa -I merged_aln/${BC}.bam -O merged_vcfs/${BC}.vcf.gz  & wait
  
		done< tmp100lines 
	i=$(echo $i+100| bc)

	echo "TOTAL PROCESSED LINES: "$i

done

rm tmp100lines 

