# MBatch R Package

This is for educational and research purposes only. 

Samples from large research projects are often processed and run in multiple batches at different times. Because the samples are processed in batches rather than all at once, the data can be vulnerable to systematic noise such as batch effects (unwanted variation between batches) and trend effects (unwanted variation over time), which can lead to misleading analysis results.

The MBatch R package is designed to help assess and correct for batch effects. It first allows the user to assess and quantify the presence of any batch effects via algorithms such as Hierarchical Clustering, Principal Component Analysis, and box plots. If significant batch effects are observed in the data, the user then has the option of selecting from a variety of correction algorithms, such as Empirical Bayes (aka Combat), ANOVA and Median Polish.

Additional information can be found at http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview

The documentation directort contains several kinds of documentation for MBatch:

 * Files that start MBatch_01 are install documentations for Linux (Debian 9.1), Windows, and OS X.
 * Files that start MBatch_02 are additional details about the test files in the package.
 * Files that start MBatch_03 are detail the file formats used by MBatch and the associated "Standardized Data" files.
 * Files that start MBatch_04 are documentation of assessment algorithms/plots.
 * Files that start MBatch_05 are documentation of correction algorithms.

Downloads and details on Standardized Data are available at http://bioinformatics.mdanderson.org/TCGA/databrowser/

See main README.MD on install instructions.

