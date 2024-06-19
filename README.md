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

- `--count_matrix COUNT_MATRIX`
 
- `--design DESIGN`

Design matrix, must have columns for ... . E.g. ./data/CTCF_demo/IN/CTCF_design_matrix_per_TF.txt.


#### **ARGUMENTS:**  

- `-h, --help`

Show this help message and exit.
  
- `--step {1,2,3,4}`

See description below.

- `--TF TF`

Should be present in the Target column in the design matrix.
  
- `--TF_label TF_LABEL`

Used for the output file prefix.
 
- `--account_covariance {TRUE,FALSE}`
  
- `--test_depletion {TRUE,FALSE}`

Whether to test for depletion instead of enrichment. Makes the values in the Cycle design matrix negative.
  
- `--step2_batches STEP2_BATCHES`

File with the name of the samples to be included in step 2.
  
- `--ltr_sample_size LTR_SAMPLE_SIZE`

How many peaks are to be included for a LTR test in step 3.
  
- `--previous_model PREVIOUS_MODEL`

Model build in step 2 to be used in step 3 or step 4.

##### **ChIP-seq arguments and files:** 

dadasd

#### **Output:**

- `--outdir OUTDIR`

Directory in which a directory named MAGIX_TF_label will be created with the output.


### **Steps**

fsfsdgfsd
sdfsdfsdf
