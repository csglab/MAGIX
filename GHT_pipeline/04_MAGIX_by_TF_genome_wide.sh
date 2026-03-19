#!/bin/bash
#SBATCH --account=def-hsn
#SBATCH --mincpus=1
#SBATCH --time=10:00:00
#SBATCH --array=
#SBATCH --mem=100000M
#SBATCH --job-name=01_MAGIX_by_TF
#SBATCH --mail-user=
#SBATCH --output=logs/%x_job_%j.out

module load StdEnv/2023 bedtools/2.31.0 r/4.3.1

if [ "${SLURM_ARRAY_TASK_ID}" = "" ]
then
    SLURM_ARRAY_TASK_ID=74
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in Login node"
    SLURM_CPUS_ON_NODE=1
else
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in working node"
fi

peaks=./ref/01_whole_genome_200_bins_wnames_wheader.bed

ref=./ref/design_matrices_by_TF_by_size.tab
wc -l ${ref}


design=./ref/by_TF/$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | \
			  awk -v FS=" " -v OFS="\t" '{print $2}' - ) 
ls -l ${design}

IN=../data/10_MAGIX_reproducibility/01_counts/v2
OUT=../data/10_MAGIX_reproducibility/02_MAGIX_by_TF_v2
experiment_id=$( basename ${design} _design_matrix.txt )_v2

OUT=${OUT}/${experiment_id}; mkdir -p ${OUT}

count_matrix=${OUT}/${experiment_id}_count_matrix_no_blacklisted.bed
DAC=./ref/ENCFF356LFX.bed 

echo ${experiment_id}; ls -l ${design}

set -o xtrace
count_files=$( tail -n +2 ${design} | \
		   awk -v FS="\t" -v OFS="\t" -v IN=${IN} '{ print IN"/"$5"_counts.txt" }' - | \
		   sed ':a;N;$!ba;s/\n/ /g' - )

target=$( tail -n +2 ${design} | \
	      awk -v FS="\t" -v OFS="\t" -v IN=${IN} '{ print $3 }' - | \
	      head -n 1 )
set +o xtrace

echo ${target}; ls -l ${peaks} ${count_files}; wc -l ${peaks} ${count_files}

paste ${peaks} ${count_files} > ${count_matrix}.tmp
head -n 1 ${count_matrix}.tmp > ${count_matrix}.tmp.header
tail -n +2 ${count_matrix}.tmp > ${count_matrix}.tmp.body

bedtools intersect -v -a ${count_matrix}.tmp.body -b ${DAC} > ${count_matrix}.tmp

cat ${count_matrix}.tmp.header ${count_matrix}.tmp > ${count_matrix}

rm ${count_matrix}.tmp.header ${count_matrix}.tmp ${count_matrix}.tmp.body

head ${count_matrix} | column -t

./src/MAGIX/MAGIX \
    --outdir ${OUT} \
    --step 2 \
    --design ${design} \
    --count_matrix ${count_matrix} \
    --TF ${target} \
    --TF_label ${experiment_id} \
    --test_depletion FALSE \
    --account_covariance TRUE
exit
