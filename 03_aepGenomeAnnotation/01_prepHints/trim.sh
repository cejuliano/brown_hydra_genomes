#! /bin/bash
#SBATCH --job-name=trim
#SBATCH -c 32
#SBATCH -t 60-0
#SBATCH --mem=36G
#SBATCH --error=trim.err
#SBATCH --output=trim.out

module load trimmomatic

prefixes=( E F3 IW M3 O W )

for arg in "${prefixes[@]}"
do
	java -jar /share/apps/Trimmomatic-0.36//trimmomatic-0.36.jar PE -threads 32 \
		reads/"$arg"_R1.fastq.gz reads/"$arg"_R2.fastq.gz \
		reads/"$arg"_R1_trim_p_fq.gz reads/"$arg"_R1_trim_up_fq.gz \
		reads/"$arg"_R2_trim_p_fq.gz reads/"$arg"_R2_trim_up_fq.gz \
		ILLUMINACLIP:./adapters.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36	
done
