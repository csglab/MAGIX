#!/bin/bash

design=../data/CTCF_demo/IN/CTCF_design_matrix_per_TF.txt
count_matrix=../data/CTCF_demo/IN/GHT_SELEX_CTCF_counts_aggregates_and_experiments_no_blacklisted_100k.bed
target=CTCF
OUT=../data/CTCF_demo/OUT

step2_model=../data/CTCF_demo/OUT/MAGIX_CTCF/target_CTCF_step_2_model.RDS
step1_selected_experiments=../data/CTCF_demo/IN/step1_selected_CTCF_experiments_per_TF.txt

echo -e "\n"; ls -l ${count_matrix} ${design}; echo -e "\n"

date; set -o xtrace
python3 ../MAGIX \
	--outdir ${OUT} \
	--step 3 \
	--count_matrix ${count_matrix} \
	--design ${design} \
	--TF ${target} \
	--TF_label ${target} \
	--previous_model ${step2_model} \
	--step2_samples ${step1_selected_experiments} \
	--lrt_sample_size 1000
set +o xtrace; date
# Demo run time < 5 min LRT of 10k peaks 

date; set -o xtrace
python3 ../MAGIX \
	--outdir ${OUT} \
	--step 3 \
	--count_matrix ${count_matrix} \
	--design ${design} \
	--TF ${target} \
	--TF_label ${target} \
	--previous_model ${step2_model} \
	--step2_samples ${step1_selected_experiments}
set +o xtrace; date

# Demo run time ~ 2 hour LRT of all 100k peaks 
exit




Mon Dec  2 02:37:26 PM EST 2024
