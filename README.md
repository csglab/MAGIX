# MAGIX: Model-based Analysis of Genomic Intervals with eXponential enrichment.

MAGIX is a statistical framework that calculates the rate of enrichment across GHT-SELEX (Genomic HT-SELEX) cycles.


#### **Requirements:** 

- Unix-compatible OS.  
- [R version 3.0.1](http://www.r-project.org/) or later.  
- R libraries: [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html), [circlize](https://jokergoo.github.io/circlize/), [ggplot2](https://www.rdocumentation.org/packages/ggplot2/versions/3.3.5), [tictoc], [Matrix], [uwot], [pheatmap], [stringr], [data.table](https://www.rdocumentation.org/packages/data.table/versions/1.14.2), [PRROC], [sads], [Rcpp], and [optparse](https://www.rdocumentation.org/packages/optparse/versions/1.6.6).
- [Python version 3.6.0](https://www.python.org/downloads/) or later.  

After cloning, you can add the line `export PATH=${DIR}/MAGIX:$PATH` to your `.bashrc` file.

#### **Usage:**  


##### **Input files:** 

- `--count_matrix COUNT_MATRIX` Count matrix each column name corresponds to the Experiment_ID column in the design file (--design).
 
 - `--design DESIGN` Design matrix, must have columns for Batch, Cycle, Experiment_ID, and Target. E.g. ./data/CTCF_demo/IN/CTCF_design_matrix_per_TF.txt.

#### **ARGUMENTS:**  

- `-h, --help` Show this help message and exit.
  
- `--step {1,2,3,4}` See description below.

- `--TF TF` Should be present in the Target column in the design matrix.
  
- `--TF_label TF_LABEL` Used for the output file prefix.
 
- `--account_covariance {TRUE,FALSE}`
  
- `--test_depletion {TRUE,FALSE}` Whether to test for depletion instead of enrichment. Makes the values in the Cycle design matrix negative.
  
- `--step2_batches STEP2_BATCHES` File with the name of the samples to be included in step 2.
  
- `--ltr_sample_size LTR_SAMPLE_SIZE` How many peaks are to be included for a LTR test in step 3.
  
- `--previous_model PREVIOUS_MODEL` Model build in step 2 to be used in step 3 or step 4.

##### **ChIP-seq arguments and files:** 

TODO

#### **Output:**

- `--outdir OUTDIR` Directory in which a directory named MAGIX_TF_label will be created with the output.


### **Steps**
Step 1: it analyzes all samples/batches for the Target TF. Creates a coefficients correlation heatmap, Q-Q plots (for all sample permutations), and estimates library size factors. 

Step 2: it limits the analysis to samples that are manually selected with the results from step 1 (`step2_batches`). Combinates the samples and fits coefficients for the Target TF.

Step 3 and 4: Performs a Likelihood-ratio test with the complete model from Step 2 and a reduced model without the Target coefficient. Step 4 performs the LRT with the entire dataset and Step 3 on a subset (size specified with `--ltr_sample_size`).





 


