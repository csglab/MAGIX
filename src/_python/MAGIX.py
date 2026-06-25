# =======================================================================
#  Authors: Aldo Hernandez-Corchado, Hamed S. Najafabadi
#
#  Copyright 2026 Aldo Hernandez-Corchado, Hamed S. Najafabadi
#
#  This file is part of MAGIX.
#
#  MAGIX is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MAGIX is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MAGIX.  If not, see <http://www.gnu.org/licenses/>.
# =======================================================================
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
