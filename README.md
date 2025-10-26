# ECV_code

The `ECV_code` repository contains all scripts and materials necessary to reproduce the results reported in our paper.  

---

## **MAIN FUNCTIONS**

The `main` directory contains the core functions used in both simulation studies and real data analyses:

- `glmcv_functions.R`: All primary R functions implementing our proposed methodology.  
- `confirm_jmle_omp_poisson_with_intercept_missing.cpp`: C++ implementation of the alternating maximization algorithm for Poisson factor models with missing data.  
- `confirm_jmle_omp_linear_with_intercept_missing.cpp`: C++ implementation of the alternating maximization algorithm for Gaussian factor models with missing data.  
- The alternating maximization algorithm for binary factor models with missing data is implemented via the **R package `mirtjml`**.

---

## **SIMULATION STUDIES**

### **SIM1 – Main Simulation Results (Examples 1–7)**

The `sim1` directory contains scripts for reproducing Examples 1–7 presented in the main text.  

Before running the code:  
- Ensure that all required R packages are installed.  
- Update the folder paths in the scripts to match the directory structure on your local machine.  
- Place the required function scripts in the specified locations to allow for proper sourcing.  
- Set the output directories to your preferred locations.  

Scripts and their corresponding results:  
- `poi_q4.R`, `output_q4.R`, `output_q4-diffpi.R`: Reproduce Table S1 (Supplementary Materials) and Example 1 in Table 1 (main text).  
- `poi_q4_othernp.R`, `output_othernp.R`: Reproduce Example 2 in Table 1.  
- `poi_q4perbp.R`, `output_perbp.R`: Reproduce Example 3 in Table 1.  
- `poi_q4weak.R`, `output_weak.R`: Reproduce Example 4 in Table 1.  
- `output_miss.R`, `output_miss.R`: Reproduce Example 5 in Table 1.  
- `poi_overdis.R`, `output_overdis.R`: Reproduce Example 6 in Table 1.  
- `poi_depen.R`, `output_depen.R`: Reproduce Example 7 in Table 1.

---

### **SIM2 – Sample Size Effects**

The `sim2` directory contains scripts for evaluating the effects of sample size (`n`) and dimensionality (`p`) on estimation performance, as presented in Figure S2 (Supplementary Materials).  

- `poi_q4_fix300p.R`, `plot_fix300p.R`: Reproduce the left panel of Figure S2 (fixed `p = 300`, varying `n`).  
- `poi_q4_fix300n.R`, `plot_fix300n.R`: Reproduce the right panel of Figure S2 (fixed `n = 300`, varying `p`).

---

### **SIM3 – Additional Simulation Studies (Examples S1–S4)**

The `sim3` directory contains scripts for supplementary simulations (Examples S1–S4, Supplementary Materials).  

- `poi_q2.R`, `output_q2.R`: Reproduce Example S1 in Table S1.  
- `poi_q4perbn.R`, `output_perbn.R`: Reproduce Example S2 in Table S1.  
- `norm_q4.R`, `output_norm.R`: Reproduce Example S3 in Table S1.  
- `bin.R`, `output_bin.R`: Reproduce Example S4 in Table S1.

---

### **SIM4 – Constraint Constant Selection**

The `sim4` directory contains scripts for analyzing the selection of the constraint constant `C`.  

- `poi_q4_diff_c.R`, `output_estc.R`, `plot_rate.R`: Reproduce the left panel of Figure S1 (results for different `C` values).  
- `estimate_c.R`: Reproduce the right panel of Figure S1 (data-driven estimation of `C`).

---

## **REAL DATA ANALYSIS**

The `real` directory contains scripts for real data analysis. To improve computational efficiency, cross-validation folds are executed separately.  

The `gene` directory contains the mouse brain dataset used in our analysis, which can also be downloaded from:  
<https://www.kaggle.com/datasets/aayush9753/singlecell-rnaseq-data-from-mouse-brain>  

Scripts:  
- `001V2real_gene8q.R`: Reproduce 5-fold cross-validation results.  
- `001V2real_gene8q-2.R`: Reproduce 10-fold cross-validation results.  
- `001V2real_gene_10miss8q-2.R`: Reproduce 10-fold cross-validation results with randomly introduced missing values.  
- `V2real_plot.R`: Generate the heatmap (Figure 1, main text) and gene importance plots (Figure S3, Supplementary Materials).

---

## **R Package: pECV**

We provide an R package, **pECV** (Penalized Entrywise Splitting Cross-Validation), to implement the proposed methodology. The package offers a user-friendly interface for applying our methods to new datasets.

### Availability

The package is available from two sources:

* **Stable Version:** The most recent stable version is available on the Comprehensive R Archive Network (CRAN). <https://cran.r-project.org/package=pECV>
* **Development Version:** The latest development version is hosted on GitHub. <https://github.com/wangATsu/ECV>

### Documentation

Full documentation, installation instructions, and usage examples are available at the GitHub repository:
https://github.com/wangATsu/ECV

---

## **COMPUTATIONAL REQUIREMENTS**

- R version 4.4.1 or higher  
- `Rcpp` package for C++ integration  
- Additional R packages as specified in each script  
- Recommended: Multi-core processor for parallel computation  


