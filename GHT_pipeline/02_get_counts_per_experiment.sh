#!/bin/bash
#SBATCH --account=def-hsn
#SBATCH --time=01:30:00
#SBATCH --mincpus=1
#SBATCH --array=20,79
#SBATCH --mem=5000M
#SBATCH --job-name=01_get_counts_per_experiments_redo2
#SBATCH --mail-user=ahcorcha@gmail.com
#SBATCH --output=logs/01_get_counts_per_experiment/%x_job_%j.out

module load StdEnv/2020 samtools/1.16.1 bedtools/2.30.0 bedops/2.4.39

if [ "${SLURM_ARRAY_TASK_ID}" = "" ]
then
    SLURM_ARRAY_TASK_ID=3
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in Login node"
    SLURM_CPUS_ON_NODE=1
else
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in working node"
fi

IN2=../data/10_MAGIX_reproducibility/filtered_bams; mkdir -p ${IN2}
OUT=../data/10_MAGIX_reproducibility/01_counts/v2; mkdir -p ${OUT}

# IN=/home/ahcorcha/projects/rrg-hsn/ahcorcha/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_SELEX/03_aligned_BAM_GHT_SELEX/
IN=../data/03_process_SELEX_redo/03_aligned_BAM_GHT_SELEX/

peaks=/home/ahcorcha/projects/rrg-hsn/ahcorcha/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_SELEX/60_whole_genome_ref/01_whole_genome_200_bins_wnames.bed

ref=./ref/Metadata_For_Aldo_And_Hamed_V22.7_AJ_AY_RR_modified_07_10_22_corrected.csv

experiment_id=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $1}' )
cycle=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $7}' )
filename=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $9}' ) 

echo ${experiment_id}; echo ${cycle}; echo ${filename}

bam=${IN}/$( basename ${filename} .fastq.gz )_filtered_sorted.bam
bam2=${IN2}/$( basename ${filename} .fastq.gz )_filtered_sorted.bam
ls -lthrg ${bam}
samtools flagstat ${bam}

count_mat=${OUT}/$( basename ${bam}  _filtered_sorted.bam )_counts.txt

echo ${count_mat}

set -o xtrace

# ### Paired bam files
samtools view -f 0x2 -b ${bam} > ${bam2}.tmp
samtools sort -n ${bam2}.tmp | samtools fixmate -m - ${bam2}.tmp2
samtools sort ${bam2}.tmp2 | samtools markdup -r - ${bam2}
samtools index ${bam2}
echo "\n\n"
samtools flagstat ${bam2}
rm ${bam2}.tmp*

### Unpaired bam files
# samtools sort -n ${bam} | samtools fixmate -m - ${bam2}.tmp2
# samtools sort ${bam2}.tmp2 | samtools markdup -r - ${bam2}
# samtools index ${bam2}
echo "\n\n"
samtools flagstat ${bam2}
rm ${bam2}.tmp*


bedtools multicov -bams ${bam2} -bed ${peaks} > ${count_mat}.tmp
set +o xtrace

head ${count_mat}.tmp; ls -l ${count_mat}.tmp

awk -v FS="\t" -v OFS="\t" '{ print $5 }' ${count_mat}.tmp | \
    sed "1s/^/${experiment_id}_cycle${cycle}\n/" - > ${count_mat}

head ${count_mat}; ls -l ${count_mat}; wc -l ${count_mat}
rm ${count_mat}.tmp
exit

5354
5354-5683


20:AATA_AffSeq_D7_ZNF384,pTH15892,ZNF384,Lysate,AffiSeqV1,D7,1,paired,ZNF384_AffSeq_Lysate_BatchAATA_Cycle1_R1.fastq.gz,ZNF384_AffSeq_Lysate_BatchAATA_Cycle1_R2.fastq.gz,Phase1
51:AATBA_AffSeq_F11_ZNF384,pTH15413,ZNF384,IVT,AffiSeqV1,F11,1,paired,ZNF384_AffSeq_IVT_BatchAATBA_Cycle1_R1.fastq.gz,ZNF384_AffSeq_IVT_BatchAATBA_Cycle1_R2.fastq.gz,Phase1
79:AATA_AffSeq_D7_ZNF384,pTH15892,ZNF384,Lysate,AffiSeqV1,D7,2,paired,ZNF384_AffSeq_Lysate_BatchAATA_Cycle2_R1.fastq.gz,ZNF384_AffSeq_Lysate_BatchAATA_Cycle2_R2.fastq.gz,Phase1
110:AATBA_AffSeq_F11_ZNF384,pTH15413,ZNF384,IVT,AffiSeqV1,F11,2,paired,ZNF384_AffSeq_IVT_BatchAATBA_Cycle2_R1.fastq.gz,ZNF384_AffSeq_IVT_BatchAATBA_Cycle2_R2.fastq.gz,Phase1
138:AATA_AffSeq_D7_ZNF384,pTH15892,ZNF384,Lysate,AffiSeqV1,D7,3,unpaired,ZNF384_AffSeq_Lysate_BatchAATA_Cycle3_R1.fastq.gz,ZNF384_AffSeq_Lysate_BatchAATA_Cycle3_R2.fastq.gz,Phase1
169



20,51,79,110,138,169



samtools view -c -F 1 ../data/03_process_SELEX_redo/03_aligned_BAM_GHT_SELEX/ZNF384_AffSeq_IVT_BatchAATBA_Cycle1_R1_filtered_sorted.bam
samtools view -c -F 1 ../data/03_process_SELEX_redo/03_aligned_BAM_GHT_SELEX/ZNF384_AffSeq_IVT_BatchAATBA_Cycle2_R1_filtered_sorted.bam
samtools view -c -F 1
