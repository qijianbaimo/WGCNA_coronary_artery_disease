## WGCNA for coronary artery disease study

* This repository is created to generate reproducible results for WGCNA network construction and hub gene identification in coronary artery disease using data from GEO ([GSE20680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20680) and [GSE20681](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20681)).

* The**download_norm.R** script is used to download and normalize data from GEO. **diff.R** script is used for further data preprocessing and **WGCNA.R** is used for downstream analysis. The scripts can be run in any operating systems (Windows, Linux, Mac, etc.) with RStudio/R and the dependent packages installed.

* The scripts may be run interactively in RStudio - an IDE for R, or run with command line **Rscript download_norm.R**, *Rscript diff.R* and *Rscript WGCNA.R*. It should be noted that the function *allowWGCNAThreads()* must be replaced by *enableWGCNAThreads()* when running the scripts in RStudio.

* It took about 10 minutes to download and normalize the raw data (*Rscript download_norm.R*), 4 minutes for data preprocessing (*Rscript diff.R*) and about 10 minutes for downstream analysis (*Rscript WGCNA.R*) on a Red Hat Linux platform with 8 cores and 16GB memory. Standard output and standard error produced by running the two scripts can be found in the **logs** directory.

* Three major input files: (1) the gene expression matrix file (**RawDataNOCtrl.txt.gz**, generated by *download_norm.R*), (2) the group information of samples used in the study (**Clinic.txt**), and (3) the probe annotation file (**GPL4133_annot.txt**), are required to run the pipeline.

* **for_WGCNA_gene.txt** is a intermediate file generated by *diff.R*, which is the input for the downstream WGCNA analysis. All the other results generated by running the pipeline can be found in the **output** directory.

* **For_gene_order.txt** file contains the order of genes in the expresion matrix for WGCNA network construction. This file can be used to generate exactly the same result as reported in the manuscript. The network maybe slightly different if not using the same order of genes according to the algorithm of WGCNA, but all the important hub genes can be identified in the modules that are most correlated with the phenotype.

* The list of MD5 checksums for all the files (input, output, and scripts, etc.) in this repository is stored in the **md5sum.txt** file. **CytoscapeInput-edges-turquoise.txt** is compressed to **CytoscapeInput-edges-turquoise.txt.gz** due to the limit of file uploading (< 100MB) in github. 


> sessionInfo()  

R version 3.5.2 (2018-12-20)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.7 (Maipo)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] readr_1.3.1           flashClust_1.01-2     WGCNA_1.68           
[4] fastcluster_1.1.25    dynamicTreeCut_1.63-1

loaded via a namespace (and not attached):
 [1] Biobase_2.42.0        bit64_0.9-7           splines_3.5.2        
 [4] foreach_1.4.4         Formula_1.2-3         assertthat_0.2.0     
 [7] stats4_3.5.2          latticeExtra_0.6-28   blob_1.2.1           
[10] fit.models_0.5-14     robustbase_0.93-3     impute_1.56.0        
[13] pillar_1.3.1          RSQLite_2.1.1         backports_1.1.3      
[16] lattice_0.20-38       glue_1.3.0            digest_0.6.23        
[19] RColorBrewer_1.1-2    checkmate_1.9.1       colorspace_1.4-1     
[22] htmltools_0.4.0       preprocessCore_1.44.0 Matrix_1.2-15        
[25] pcaPP_1.9-73          pkgconfig_2.0.2       purrr_0.3.3          
[28] GO.db_3.7.0           mvtnorm_1.0-12        scales_1.0.0         
[31] htmlTable_1.13.1      tibble_2.1.3          IRanges_2.16.0       
[34] ggplot2_3.2.1         nnet_7.3-12           BiocGenerics_0.28.0  
[37] lazyeval_0.2.1        survival_2.43-3       magrittr_1.5         
[40] crayon_1.3.4          memoise_1.1.0         doParallel_1.0.14    
[43] MASS_7.3-51.1         foreign_0.8-71        tools_3.5.2          
[46] data.table_1.12.0     hms_0.4.2             matrixStats_0.54.0   
[49] stringr_1.4.0         S4Vectors_0.20.1      munsell_0.5.0        
[52] cluster_2.0.7-1       AnnotationDbi_1.44.0  compiler_3.5.2       
[55] rlang_0.4.4           grid_3.5.2            iterators_1.0.10     
[58] rstudioapi_0.10       htmlwidgets_1.3       robust_0.4-18.1      
[61] base64enc_0.1-3       gtable_0.2.0          codetools_0.2-15     
[64] DBI_1.0.0             rrcov_1.4-7           R6_2.4.0             
[67] gridExtra_2.3         knitr_1.27            dplyr_0.8.3          
[70] bit_1.1-14            Hmisc_4.2-0           stringi_1.3.1        
[73] parallel_3.5.2        Rcpp_1.0.1            vctrs_0.2.2          
[76] rpart_4.1-13          acepack_1.4.1         DEoptimR_1.0-8       
[79] tidyselect_0.2.5      xfun_0.8 
