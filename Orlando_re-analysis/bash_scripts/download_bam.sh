#!/bin/bash
#SBATCH --job-name=rna_seq_bam_dowload
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=8000M
#SBATCH -p long
#SBATCH --output=rna_seq_dowload.out
#SBATCH --error=rna_seq_dowload.error
#SBATCH --mail-user=m.corteslopez@imb-mainz.de
#SBATCH --mail-type=ALL

module load sratoolkit
while read line;
do
	sam-dump ${line} | samtools view -bS - > ${line}.bam
done< ../SRA_IDs.list



