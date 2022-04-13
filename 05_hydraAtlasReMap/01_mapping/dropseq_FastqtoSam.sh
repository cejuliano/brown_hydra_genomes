#!/bin/bash -l
#SBATCH --job-name=ds
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --array=1-11
#SBATCH --error=logs/fstqsam.%A_%a.err
#SBATCH --output=logs/fstqsam.%A_%a.out
#SBATCH -t 2-00:00:00

# ======================================== #
### DEFINE PATHS NEEDED FOR THE PIPELINE ###
# ======================================== #

# Path to .jar files for Picard tools
#picard_path="/seq/software/picard/current/bin/"

export PATH=/group/julianolab/analyses/dropseq/jre1.8.0_121/bin:$PATH

module load bowtie2
module load samtools

picard_path="/group/julianolab/analyses/dropseq/Drop-seq_tools-2.4.0/3rdParty/picard"

# Path to Hydra genome metadata (bowtie indexes, fasta files,
# .gtf, .refFlat, etc.)
# NOTE: paths include prefix of these files as well as directory.

# Path to dropseq tools

dropseq_tools="/group/julianolab/analyses/dropseq/Drop-seq_tools-2.4.0"

# Path to files to work on (configure these)
fastq_dir="/group/julianolab/data/hydra_ds/" #fastq directory
main_dir="/group/julianolab/analyses/dropseq/HVAEP1_transcriptome_final/" 
temp_dir="/group/julianolab/analyses/dropseq/HVAEP1_transcriptome_final/temp/"
file_list="/group/julianolab/analyses/dropseq/HVAEP1_transcriptome_final/file_list_all.txt"

# Paths that are automatically determined
#file_prefix=`sed -n "$SGE_TASK_ID"p "${file_list}" | cut -f 1`
file_prefix=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1`
work_dir="${main_dir}${file_prefix}/"
temp_work_dir="${temp_dir}${file_prefix}/"
fastq_path="${fastq_dir}/${file_prefix}"
log_dir="${work_dir}/logs/"

# What sections of the pipeline do you want to run?
# 0 = OFF
# 1 = ON
run_1=1  #  1: Convert FASTQ to SAM

# =================================== #
### DEFINE FUNCTIONS USED IN SCRIPT ###
# =================================== #

# Echo command and then execute it.
function v_exe
{
    echo "$1"
    eval "$1" || error_exit "Cannot execute command: $1"
}
# Echo a message to both stdout and stderr to mark parts of pipeline
function echo_both
{
	echo "$1"
	echo "$1" > /dev/stderr
}
# Function to preserve the logs after each step
function save_logs
{
	v_exe "cp ${main_dir}/logs/*.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${log_dir}/${file_prefix}.${1}.out"
	v_exe "cp ${main_dir}/logs/*.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${log_dir}/${file_prefix}.${1}.err"
}


# ================================= #
### PIPELINE ACTION IS BELOW HERE ###
# ================================= #

# Get the work & log directories ready
mkdir -p ${log_dir}
#rm ${log_dir}/*
v_exe "mkdir -p ${temp_dir}"
v_exe "mkdir -p ${temp_work_dir}"
v_exe "mkdir -p ${work_dir}"

# First, we have to convert the raw FASTQ files to SAM format for the Dropseq tools
# to work on.
if [ "${run_1}" -eq 1 ]; then
	echo_both ""
	echo_both "----- 1: FASTQ to SAM -----"
	v_exe "java -Xmx4g -jar ${picard_path}/picard.jar FastqToSam FASTQ=${fastq_path}_R1_001.fastq.gz FASTQ2=${fastq_path}_R2_001.fastq.gz SAMPLE_NAME='${file_prefix}' OUTPUT=${temp_work_dir}/${file_prefix}.bam"
	save_logs "01"
fi
