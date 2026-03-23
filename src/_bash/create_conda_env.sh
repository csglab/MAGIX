#!/bin/bash

conda create --name MAGIX_env r-base r-essentials
conda activate MAGIX_env

which R Rscript

Rscript ./src/_R/install_R_lbraries.R
