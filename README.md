# MAGIX: Model-based Analysis of Genomic Intervals with eXponential enrichment.

MAGIX is a generative model that explicitly connects the enrichment of TF-bound genomic intervals to the fragment counts observed across GHT-SELEX cycles. MAGIX, models how TF-bound intervals progressively occupy a higher proportion of selected fragments pool in each cycle relative to genomic background. These fragment proportions, in turn, are treated as latent variables in the model that, together with a sample-specific library size factor, determine the number of observed reads through a Poisson process.


## **Requirements** 

- Unix-compatible OS.  
- [R version 3.0.1](http://www.r-project.org/) or later.  
- R libraries: [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html), [circlize](https://jokergoo.github.io/circlize/), [ggplot2](https://www.rdocumentation.org/packages/ggplot2/versions/3.3.5), [tictoc](https://cran.r-project.org/web/packages/tictoc/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [uwot](https://cran.r-project.org/web/packages/uwot/index.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html), [data.table](https://www.rdocumentation.org/packages/data.table/versions/1.14.2), [PRROC](https://cran.r-project.org/web/packages/PRROC/index.html), [sads](https://cran.r-project.org/web/packages/sads/index.html), [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html), and [optparse](https://www.rdocumentation.org/packages/optparse/versions/1.6.6).
- [Python version 3.6.0](https://www.python.org/downloads/) or later.  


## **Installation** 

```bash
git clone https://github.com/csglab/MAGIX.git
```

After cloning, you can add the line `export PATH=${cloning_directory}/MAGIX:$PATH` to your `.bashrc` file.

## **Usage**  

**Step 1**: Analyzes every sample for the Target TF separately. Creates a coefficients correlation heatmap, Q-Q plots (for all sample permutations), and estimates library size factors. 

**Step 2**: Combines the samples and fits coefficients for the target TF itself. You can limit the analysis to samples that are selected with the results from step 1 (e.g. `*_step_1_coefs_heatmap.pdf`) using `--step2_batches` (default = `all`). 

**Step 3**: Performs a Likelihood-ratio test with the complete model from `Step 2` and a reduced model without the target TF coefficient. It performs the LRT with the entire dataset by default, can be limitated to a smaller sample with `--ltr_sample_size`.


## **Demo**
Demo scripts and corresponding dataset can be found in `./demo_scripts` and `./data/CTCF_demo/IN`. Demo run time for `Step 1` and `2` is ~ 5 min each and ~2 hours for `Step 3`. 

For example:
```bash
MAGIX --outdir ${out_directory} \
        --step 1 \
        --count_matrix ${count_matrix} \
        --design ${design_file} \
        --TF CTCF \
        --TF_label CTCF_FL \
        --account_covariance TRUE \
        --test_depletion FALSE
```

You can download the demo output from: https://usegalaxy.org/u/ahcorcha/h/magixdemooutput.

## **Input files** 

- `--count_matrix COUNT_MATRIX` Count matrix each column name corresponds to the Experiment_ID column in the design file (--design).
 - `--design DESIGN` Design matrix, must have columns for Batch, Cycle, Experiment_ID, and Target. E.g. `./data/CTCF_demo/IN/CTCF_design_matrix_per_TF.txt`.
 - `--step2_batches` For step 2 and 3, sample names to include in the analysis (default = `all`). 


## **Citation**

Jolma, A., Hernandez-Corchado, A., Yang, A. W. H., Fathi, A., Laverty, K. U., Brechalov, A., Razavi, R., Albu, M., Zheng, H., Kulakovskiy, I. V., Najafabadi, H. S., & Hughes, T. R. (2024). GHT-SELEX demonstrates unexpectedly high intrinsic sequence specificity and complex DNA binding of many human transcription factors. BioRxiv, 2024.11.11.618478. https://doi.org/10.1101/2024.11.11.618478





Heatmaps are generated with the ComplexHeatmap and circlize packages. If you use them in published research, please cite:

- Gu, Z. Complex Heatmap Visualization. iMeta 2022.
or
- Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional
    genomic data. Bioinformatics 2016.
and
- Gu, Z. circlize implements and enhances circular visualization
  in R. Bioinformatics 2014.





 


