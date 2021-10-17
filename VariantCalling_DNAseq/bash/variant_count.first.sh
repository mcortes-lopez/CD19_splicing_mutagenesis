#!/bin/bash
#SBATCH --job-name=count_vcf_first
#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --mem=4000M
#SBATCH -p short
#SBATCH --output=count_vcf_first.out
#SBATCH --error=count_vcf_first.error
#SBATCH --mail-user=m.corteslopez@imb-mainz.de
#SBATCH --mail-type=ALL


find merged_vcfs/ -name "*vcf.gz" > TMP_FILES

totallines=$(ls merged_vcfs/| grep "gz$" | wc -l)
i=1
while [ $i -le $totallines ]
	do
	awk -v st="$i" 'NR>=st && NR< st+100{print $0}' TMP_FILES > tmp100lines
	while read line; 
		do
		BC=$line
		variants=$(zcat $line | grep -c -v "#")
		echo -e "${BC}\t${variants}" | sed 's/merged_vcfs\///g;s/.vcf.gz//g' >> variant_count.txt &
	done< tmp100lines 
	i=$(echo $i+100| bc)
done

rm tmp100lines 
rm TMP_FILES

