#!/bin/bash
#SBATCH --account=def-hsn
#SBATCH --mincpus=1
#SBATCH --time=05:30:00
#SBATCH --array=5,7
#SBATCH --mem=20000M
#SBATCH --job-name=MAGIX_peaks_version_2
#SBATCH --mail-user=ahcorcha@gmail.com
#SBATCH --output=logs/MAGIX_peaks_version_2/%x_job_%j.out

if [ "${SLURM_ARRAY_TASK_ID}" = "" ]
then
    SLURM_ARRAY_TASK_ID=86
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in Login node"
    SLURM_CPUS_ON_NODE=1
else
    echo "SLURM_ARRAY_TASK_ID: "${SLURM_ARRAY_TASK_ID}
    echo "Running in working node"
fi

round() {
    local df=${2:-3}
    printf '%.*f\n' "$df" "$(bc -l <<< "a=$1; if(a>0) a+=5/10^($df+1) else if (a<0) a-=5/10^($df+1); scale=$df; a/1")"
}

module load StdEnv/2023 bedtools/2.31.0 r/4.3.1 samtools/1.20

ref=./ref/design_matrices_by_TF_by_size.tab
wc -l ${ref}

design=./ref/by_TF/$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ref} | \
			  awk -v FS=" " -v OFS="\t" '{print $2}' - ) 
ls -l ${design}

IN=../data/10_MAGIX_reproducibility/filtered_bams
OUT=../data/10_MAGIX_reproducibility/MAGIX_peaks_version_2

experiment_id=$( basename ${design} _design_matrix.txt )

OUT=${OUT}/${experiment_id}; mkdir -p ${OUT}

count_matrix=${OUT}/${experiment_id}_count_matrix_no_blacklisted.bed

echo ${experiment_id}; ls -l ${design}

set -o xtrace
bam_files=$( tail -n +2 ${design} | \
		 sed '/Aggreg/d' - | \
		   awk -v FS="\t" -v OFS="\t" -v IN=${IN} '{ print IN"/"$5"*.bam" }' - | \
		   sed ':a;N;$!ba;s/\n/ /g' - )

target=$( tail -n +2 ${design} | \
	      awk -v FS="\t" -v OFS="\t" -v IN=${IN} '{ print $3 }' - | \
	      head -n 1 )

ls -lthrg ${bam_files}
ls -1 ${bam_files} | wc -l

prefix=${OUT}/${experiment_id}

########################################################################################## ##
samtools merge -f -o ${prefix}_merged.bam ${bam_files}
samtools index ${prefix}_merged.bam
samtools sort -T ${prefix}.tmp. -n ${prefix}_merged.bam > ${prefix}_merged_sorted_by_name.bam

properly_paired=$( samtools view -c -f 1 ${prefix}_merged_sorted_by_name.bam )


echo -e "\n\n"
samtools view -c -f 1 ${prefix}_merged_sorted_by_name.bam
samtools view -c -F 1 ${prefix}_merged_sorted_by_name.bam
echo -e "\n\n"

ls -lthrg ${prefix}_merged_sorted_by_name.bam

samtools flagstat ${prefix}_merged_sorted_by_name.bam


########################################################################################## ##
echo -e "\nUnpaired reads\n"
samtools view -b -F 1 ${prefix}_merged_sorted_by_name.bam | \
    bedtools bamtobed -i - | \
    awk -v FS="\t" -v OFS="\t" '{ print $1,$2,$3,1}' - | \
    sed '/chrM/d' - | \
    sort -k 1,1 -k2,2n - > ${prefix}_fragments.bed.tmp
head ${prefix}_fragments.bed.tmp


echo -e "\nPaired reads\n"
samtools view -b -f 1 ${prefix}_merged_sorted_by_name.bam | \
    bedtools bamtobed -bedpe -i - | \
    awk -v FS="\t" -v OFS="\t" '{ print $1,$2,$6,1}' - | \
    sed '/chrM/d' - | \
    sort -k 1,1 -k2,2n - >> ${prefix}_fragments.bed.tmp

sort -k 1,1 -k2,2n ${prefix}_fragments.bed.tmp > ${prefix}_fragments.bed

head ${prefix}_fragments.bed


bedtools merge -i ${prefix}_fragments.bed > ${prefix}_fragments_merged.bed

bedtools coverage -a ${prefix}_fragments_merged.bed -b ${prefix}_fragments.bed -d > \
         ${prefix}_coverage.bed

lib_size=$( samtools view -c ${prefix}_merged.bam )

cutoff=$( round "(${lib_size})/2000000" 1 | bc | awk '{print int($1)}' )

if [ $(( ${cutoff} )) -le 4  ];then
    cutoff=4
fi


echo ${experiment_id}" Cutoff: "${cutoff}

head -n 3 ${prefix}_coverage.bed

awk -v FS="\t" -v OFS="\t" '{ print $1,$2+$4-1,$2+$4,$5}' ${prefix}_coverage.bed > \
    ${prefix}_coverage_formated.bed

head -n 3 ${prefix}_coverage_formated.bed

bedtools merge -c 4 -o max -i ${prefix}_coverage_formated.bed | \
    awk -v FS="\t" -v OFS="\t" -v cutoff=${cutoff} \
	'{ if ($4 >= cutoff) print $1,$2,$3,$4}' - > \
	${prefix}_coverage_formated_filtered.bed

head -n 3 ${prefix}_coverage_formated_filtered.bed
wc -l ${prefix}_coverage_formated_filtered.bed


bedtools intersect -wo -a ${prefix}_coverage_formated.bed -b ${prefix}_coverage_formated_filtered.bed | \
    awk -v FS="\t" -v OFS="\t" '{ print $1,$2,$3, $5":"$6"-"$7 ,$4}' - > \
	${prefix}_coverage_filtered_merged_both.bed

head -n 3 ${prefix}_coverage_filtered_merged_both.bed
wc -l ${prefix}_coverage_filtered_merged_both.bed

Rscript ./src/get_summits.R --in_bed ${prefix}_coverage_filtered_merged_both.bed --out_bed ${prefix}_summit.bed

awk -v FS="\t" -v OFS="\t" '{ print $1,$2-100,$2+100 }' ${prefix}_summit.bed | \
    awk -v FS="\t" -v OFS="\t" '{ print $1,$2,$3,$1":"$2"-"$3 }' -  > ${prefix}_peaks.bed

head ${prefix}_peaks.bed
wc -l ${prefix}_peaks.bed

set -o xtrace

##############################################################################################
##############################################################################################
mat=${prefix}_count_matrix.tab

list_of_BAMS=$( tail -n +2 ${design} | \
		    awk -v FS="\t" -v OFS="\t" -v BAMS=${IN} \
			'{ print BAMS"/"$5"*.bam "}' - | \
		    tr -d '\n' )

tail -n +2 ${design} | \
    awk -v FS="\t" -v OFS="\t" -v BAMS=${IN} '{ print $1}' - | \
    sed '1s|^|chr\tstart\tstop\tname\t|' - | \
    tr '\n' '\t' | \
    sed '1s|\t$|\n|' - > ${mat}

bedtools multicov -bams ${list_of_BAMS} -bed ${prefix}_peaks.bed >> ${mat}

step2_model=../data/10_MAGIX_reproducibility/02_MAGIX_by_TF_v2/${experiment_id}_v2/MAGIX_${experiment_id}_v2/target_${experiment_id}_v2_step_2_model.RDS

ls -l ${step2_model}

#############################################################################################
echo -e "\n\n################################### LRT \n"
set -o xtrace
./src/MAGIX/MAGIX \
    --outdir ${OUT} \
    --step 3 \
    --count_matrix ${mat} \
    --design ${design} \
    --TF ${experiment_id} \
    --TF_label ${experiment_id} \
    --account_covariance TRUE \
    --test_depletion FALSE \
    --previous_model ${step2_model}
set +o xtrace

head ${OUT}/MAGIX_${experiment_id}/target_${experiment_id}_step_3_LTR_results.csv

set -o xtrace
echo -e "\n\n"
tail -n +2 ${OUT}/MAGIX_${experiment_id}/target_${experiment_id}_step_3_LTR_results.csv | \
    sed 's|"||g' - | \
    sed 's|-|,|' - | \
    sed 's|:|,|' - | \
    awk -v FS="," -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3,$4,$5,$6,$7,$8,$9}' - | \
    sort -k5,5gr - > \
	 ${OUT}/MAGIX_${experiment_id}/target_${experiment_id}_step_3_LTR_results.bed

head ${OUT}/MAGIX_${experiment_id}/target_${experiment_id}_step_3_LTR_results.bed

echo "chr,start,stop,name,coefficient.br,coefficient.ar,full_LL,reduced_LL,pvalue,fdr" > \
     ${OUT}/MAGIX_${experiment_id}/HEADER_for_step_3_LTR_results.bed


# rm ${prefix}_merged.bam ${prefix}_merged_sorted.bam
exit















1,3-5,7-9,11-206


2h
1,2,3,4,5,6,7,8,10,11,15,19,23,24,26,28,31,35,49,50,55,61,74,75,117

30m
9,12-14,16-18,20-22,25,27,29,30,32-34,36-48,51-54,56-60,62-73,76-116,118-206




# samtools sort ${prefix}_merged.bam > ${prefix}_merged_sorted_by_coord.bam
# samtools index ${prefix}_merged_sorted_by_coord.bam
DAC=./ref/ENCFF356LFX.bed 



-rw-r----- 1 rrg-hsn  286K Aug 15 23:54 ./SNAI1/MAGIX_SNAI1/target_SNAI1_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn   47K Aug 15 23:55 ./ZNF322/MAGIX_ZNF322/target_ZNF322_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn   47K Aug 15 23:55 ./ZNF260/MAGIX_ZNF260/target_ZNF260_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn   43K Aug 15 23:56 ./ZNF596/MAGIX_ZNF596/target_ZNF596_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  384K Aug 15 23:57 ./ZNF436/MAGIX_ZNF436/target_ZNF436_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  354K Aug 15 23:58 ./ZFP3/MAGIX_ZFP3/target_ZFP3_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  1.9M Aug 15 23:58 ./ZNF18/MAGIX_ZNF18/target_ZNF18_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  1.8M Aug 15 23:58 ./ZIM3/MAGIX_ZIM3/target_ZIM3_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn   89K Aug 15 23:58 ./ZNF264/MAGIX_ZNF264/target_ZNF264_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  1.2M Aug 15 23:59 ./ZNF134/MAGIX_ZNF134/target_ZNF134_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  373K Aug 16 00:00 ./ZNF770/MAGIX_ZNF770/target_ZNF770_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  359K Aug 16 00:03 ./ZNF16/MAGIX_ZNF16/target_ZNF16_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  205K Aug 16 00:03 ./GLI4/MAGIX_GLI4/target_GLI4_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  731K Aug 16 00:03 ./ZNF250/MAGIX_ZNF250/target_ZNF250_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn  1.5M Aug 16 00:09 ./ZNF35/MAGIX_ZNF35/target_ZNF35_step_3_LTR_results.bed
-rw-r----- 1 rrg-hsn   18M Aug 16 01:11 ./ZNF384/MAGIX_ZNF384/target_ZNF384_step_3_LTR_results.bed


301 ./ZNF596/MAGIX_ZNF596/target_ZNF596_step_3_LTR_results.bed
319 ./ZNF260/MAGIX_ZNF260/target_ZNF260_step_3_LTR_results.bed
374 ./ZNF322/MAGIX_ZNF322/target_ZNF322_step_3_LTR_results.bed
651 ./ZNF264/MAGIX_ZNF264/target_ZNF264_step_3_LTR_results.bed



## ATAA and ATAB
#SBATCH --array=74,102,103,108,141,162,174,178,179,180,188,189,192,195,203,206

5,7,86
5  - FLI1_design_matrix.txt
86 - MSANTD1_design_matrix.txt
7  - MYPOP_design_matrix.txt





# if [ ${properly_paired} = 0  ];then

#     echo -e "\nUnpaired reads\n"
#     bedtools bamtobed -i ${prefix}_merged_sorted_by_name.bam | \
# 	awk -v FS="\t" -v OFS="\t" '{ print $1,$2,$3,1}' - | \
# 	sed '/chrM/d' - | \
# 	sort -k 1,1 -k2,2n - > ${prefix}_fragments.bed
#     head ${prefix}_fragments.bed
# fi

# if [ ${properly_paired} != 0  ];then

#     echo -e "\nPaired reads\n"
#     bedtools bamtobed -bedpe -i ${prefix}_merged_sorted_by_name.bam | \
# 	awk -v FS="\t" -v OFS="\t" '{ print $1,$2,$6,1}' - | \
# 	sed '/chrM/d' - | \
# 	sort -k 1,1 -k2,2n - > ${prefix}_fragments.bed
#     head ${prefix}_fragments.bed
# fi
