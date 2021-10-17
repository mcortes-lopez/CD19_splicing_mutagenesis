#!/bin/bash
#SBATCH --job-name=rna_seq_dowload
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=4000M
#SBATCH -p long
#SBATCH --output=rna_seq_dowload.out
#SBATCH --error=rna_seq_dowload.error
#SBATCH --mail-user=m.corteslopez@imb-mainz.de
#SBATCH --mail-type=ALL

module load sratoolkit
while read line;
do
	fastq-dump.2.8.1-3 --outdir ../fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --accept-hard-clip ${line}
done< ../SRA_IDs.list



