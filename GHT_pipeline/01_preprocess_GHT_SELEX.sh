#!/bin/bash
#SBATCH --account=def-hsn
#SBATCH --time=00:30:00
#SBATCH --array=
#SBATCH --mincpus=1
#SBATCH --mem=4000M
#SBATCH --job-name=01_align_GHT_SELEX
#SBATCH --mail-user=
#SBATCH --output=logs/%x_job_%j.out

if [ "${SLURM_ARRAY_TASK_ID}" = "" ]
then
    SLURM_ARRAY_TASK_ID=178 # 136 # 2-5684
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in Login node"
    SLURM_CPUS_ON_NODE=1
else
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in working node"
fi

module load StdEnv/2020 trimmomatic/0.39 bowtie2/2.4.4 fastqc/0.11.9 samtools/1.16.1


HG=/home/ahcorcha/projects/rrg-hsn/ahcorcha/ahcorcha/resources/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
ref=./ref/Metadata_For_Aldo_And_Hamed_V22.7_AJ_AY_RR_modified_07_10_22_corrected.csv
IN=../data/03_process_SELEX/01_FQ_original_ajolma_SELEX/RawData

exp_id=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $1}' )
cycle=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $7}' )
read_type=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $8}' )
sub_dir=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $11}' )

R1=${IN}/${sub_dir}/$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $9}' )
R2=${IN}/${sub_dir}/$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | awk -v FS="," '{print $10}' )

OUT_TRIMM=../data/03_process_SELEX_redo/02_FQ_trimmed_GHT_SELEX
OUT_BAM=../data/03_process_SELEX_redo/03_aligned_BAM_GHT_SELEX

mkdir -p ${OUT_BAM} ${OUT_TRIMM}
basename=$( basename ${R1} .fastq.gz )


R1_trim=${OUT_TRIMM}/${basename}_trimmed_R1.fastq.gz
R1_U=${OUT_TRIMM}/${basename}_unpaired_R1.fastq.gz

R2_trim=${OUT_TRIMM}/${basename}_trimmed_R2.fastq.gz
R2_U=${OUT_TRIMM}/${basename}_unpaired_R2.fastq.gz

stats=${OUT_TRIMM}/${basename}_statsSummary.txt

BAM_unfilt=${OUT_BAM}/${basename}_unfiltered.bam
BAM_filt=${OUT_BAM}/${basename}_filtered.bam
BAM_filt_sort=${OUT_BAM}/${basename}_filtered_sorted.bam


echo "Experiment ID: "${exp_id}
echo "        Cycle: "${cycle}
echo "    Read type: "${read_type}
echo " R1 file name: "${R1}
echo " R2 file name: "${R2}
echo "Sub directory: "${sub_dir}
echo "     Basename: "${basename}
echo -e "\n"

if [ "${read_type}" = "paired" ]
then
    ls -l ${R1} ${R2}; echo -e "\n"

    set -o xtrace
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
	 PE -threads ${SLURM_CPUS_ON_NODE} \
	 -summary ${stats} \
	 ${R1} ${R2}  \
	 ${R1_trim} \
	 ${R1_U} \
	 ${R2_trim} \
	 ${R2_U} \
	 ILLUMINACLIP:ref/adapters_CustomAffiseqSE.fa:2:5:5 LEADING:3 TRAILING:3 MINLEN:25
    set +o xtrace
    
    echo -e "\n"; ls -lShg ${R1_trim} ${R1_U} ${R2_trim} ${R2_U}; echo -e "\n"
    echo "######## Trimming stats"; cat ${stats}; echo -e "\n"
    
    ## ReadMapping
    set -o xtrace
    bowtie2 --very-sensitive --no-unal -x ${HG} \
	    --threads ${SLURM_CPUS_ON_NODE} \
	    -1 ${R1_trim} -2 ${R2_trim} \
	| samtools view --threads ${SLURM_CPUS_ON_NODE} -bh -o ${BAM_unfilt} -
    set +o xtrace
    
    rm ${R1_U} ${R2_U} ${R1_trim} ${R2_trim}
fi

if [ "${read_type}" = "unpaired" ]
  then
    ls -lthg ${R1}; echo -e "\n"

    set -o xtrace
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
	 SE -threads ${SLURM_CPUS_ON_NODE} \
	 -summary ${stats} \
	 ${R1} ${R1_trim} \
	 ILLUMINACLIP:./ref/adapters_CustomAffiseqSE.fa:2:5:5 LEADING:3 TRAILING:3 MINLEN:25
    set +o xtrace
    
    echo -e "\n"; ls -lShg ${R1_trim}; echo -e "\n"
    echo "######## Trimming stats"; cat ${stats}; echo -e "\n"
    
    ## ReadMapping
    set -o xtrace
    bowtie2 --very-sensitive --no-unal -x ${HG} \
	    --threads ${SLURM_CPUS_ON_NODE} \
	    -U ${R1_trim} \
	| samtools view --threads ${SLURM_CPUS_ON_NODE} -bh -o ${BAM_unfilt} -
    set +o xtrace
    
    rm ${R1_trim}
fi

## Filtering
set -o xtrace
samtools view --threads ${SLURM_CPUS_ON_NODE} -bhq 30 -F 1548 ${BAM_unfilt} > ${BAM_filt}
## No MAPQ cutoff
# samtools view --threads ${SLURM_CPUS_ON_NODE} -bh -F 1548 ${BAM_unfilt} > ${BAM_filt}

samtools sort --threads ${SLURM_CPUS_ON_NODE} ${BAM_filt} -o ${BAM_filt_sort}
samtools index -@ ${SLURM_CPUS_ON_NODE} ${BAM_filt_sort}



ls -lthgS ${BAM_filt} ${BAM_unfilt} ${BAM_filt_sort}
rm ${BAM_filt} ${BAM_unfilt}

exit
