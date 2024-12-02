#!/bin/bash

design=../data/CTCF_demo/IN/CTCF_design_matrix_per_TF.txt
count_matrix=../data/CTCF_demo/IN/GHT_SELEX_CTCF_counts_aggregates_and_experiments_no_blacklisted_100k.bed
target=CTCF
OUT=../data/CTCF_demo/OUT

step1_selected_experiments=../data/CTCF_demo/IN/step1_selected_CTCF_experiments_per_TF.txt

echo -e "\n"; ls -l ${count_matrix} ${design}; echo -e "\n"

date; set -o xtrace
python3 ../MAGIX \
	--outdir ${OUT} \
	--step 2 \
	--count_matrix ${count_matrix} \
	--design ${design} \
	--TF ${target} \
	--TF_label ${target} \
	--account_covariance TRUE \
	--test_depletion FALSE \
	--step2_samples ${step1_selected_experiments}
set +o xtrace; date
exit

# Demo run time < 5 min
