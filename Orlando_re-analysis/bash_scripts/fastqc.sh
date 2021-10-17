#!/bin/bash
#SBATCH --job-name=quality_check
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=4000M
#SBATCH -p short
#SBATCH --output=quality_check.out
#SBATCH --error=quality_check.error
#SBATCH --mail-user=m.corteslopez@imb-mainz.de
#SBATCH --mail-type=ALL

module load fastqc
find ./ -name *.fastq.gz | while read line;
do
	fastqc ${line} -o fastqc/ -t 4
done 



