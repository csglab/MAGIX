#!/bin/bash
#SBATCH --account=def-hsn
#SBATCH --time=00:30:00
#SBATCH --mincpus=1
#SBATCH --array=1,20,39
#SBATCH --mem=10000M
#SBATCH --job-name=02_get_counts_per_experiment_aggregates_v2
#SBATCH --mail-user=ahcorcha@gmail.com
#SBATCH --output=logs/02_get_counts_per_experiment_aggregates/%x_job_%j.out

if [ "${SLURM_ARRAY_TASK_ID}" = "" ]
then
    SLURM_ARRAY_TASK_ID=1
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in Login node"
else
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in working node"
fi

OUT=../data/10_MAGIX_reproducibility/01_counts/v2/; mkdir -p ${OUT}

ref=./ref/files_per_batch.txt
## ls -1 ./ref/files_per_batch/*_list.txt > ${ref}; exit

# ref=./ref/bam_files_per_batch.txt
## ls -1 ./ref/bam_files_per_batch/*txt > ${ref}; exit

batch_files=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} ) 
wc -l ${batch_files}

name=$( basename ${batch_files} _list.txt )
echo ${name}
cat ${batch_files}

echo -e "\n\n"
BAM_DIR=../data/10_MAGIX_reproducibility/filtered_bams/
BAM=${BAM_DIR}/$( basename ${batch_files} _list.txt )_aggregate.bam
echo ${BAM}

# module load StdEnv/2020 samtools/1.16.1
# set -o xtrace
# samtools merge -f --threads 1 -o ${BAM} -b ${batch_files}
# samtools index ${BAM}
# set +o xtrace
# exit

echo -e "\n\n"
module load StdEnv/2023 r/4.3.1
set -o xtrace
Rscript ./src/sum_batch_counts.R \
	--file_with_filenames ${batch_files} \
	--header ${name} \
	--outfile ${OUT}/${name}_counts.txt
exit

## cycle_1_batches_YWK_D_list.txt
## cycle_3_batches_YWK_D_counts.txt
## title: cycle_1_batches_AATA
# YWL_B YWC_C AATA
1,20,39,2,21,40
