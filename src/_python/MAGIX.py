import csv   
import os
import subprocess

def check_source_files( args, log ):
    pass

def run_cmd( cmd, log ):
    print( cmd )
    subprocess.call( cmd, shell = True )
    log.write("%s\n" % ( cmd ))

def run_MAGIX( args, log ):

    cmdline = f"""Rscript {args.}/src/_R/MAGIX.R --script_path {args.} --outdir {args.OUTDIR} --step {args.STEP} --count_matrix {args.COUNT_MATRIX} --design {args.DESIGN} --TF {args.TF} --TF_label {args.TF_LABEL} --step2_batches {args.STEP2_SAMPLES} --account_covariance {args.ACCOUNT_COVARIANCE} --ltr_sample_size {args.LTR_SAMPLE_SIZE} --test_depletion {args.TEST_DEPLETION} --use_chip {args.USE_CHIP} --chip_peaks {args.CHIP_PEAKS} --only_compare_to_chip {args.ONLY_COMPARE_TO_CHIP} --total_ChIP {args.TOTAL_CHIP} """

    utils.run_cmd(cmdline, log)
