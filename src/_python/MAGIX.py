import csv   
import os
import subprocess

def check_source_files( args ):
    pass

def run_cmd( cmd ):
    print( cmd )
    subprocess.call( cmd, shell = True )
    # log.write("%s\n" % ( cmd ))

def run_MAGIX( args ):
    
    cmdline = f"""Rscript {args.script_path}/src/_R/MAGIX.R --script_path {args.script_path} --outdir {args.OUTDIR} --step {args.STEP} --count_matrix {args.COUNT_MATRIX} --design {args.DESIGN} --TF {args.TF} --TF_label {args.TF_LABEL} --step2_batches {args.STEP2_SAMPLES} --account_covariance {args.ACCOUNT_COVARIANCE} --lrt_sample_size {args.LRT_SAMPLE_SIZE} --test_depletion {args.TEST_DEPLETION} --use_chip {args.USE_CHIP} --chip_peaks {args.CHIP_PEAKS} --only_compare_to_chip {args.ONLY_COMPARE_TO_CHIP} --total_ChIP {args.TOTAL_CHIP} --previous_model {args.PREVIOUS_MODEL} --remove_target_from_aggregate {args.RM_FROM_AGG} """

    run_cmd( cmdline )
