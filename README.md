# MBatch and MBatchUtils R Packages and MBatch Stand-Alone Docker Image with RStudio

This is for educational and research purposes only. 

Samples from large research projects are often processed and run in multiple batches at different times. Because the samples are processed in batches rather than all at once, the data can be vulnerable to systematic noise such as batch effects (unwanted variation between batches) and trend effects (unwanted variation over time), which can lead to misleading analysis results.

The MBatch R package is designed to help assess and correct for batch effects. It first allows the user to assess and quantify the presence of any batch effects via algorithms such as Hierarchical Clustering, Principal Component Analysis, and box plots. If significant batch effects are observed in the data, the user then has the option of selecting from a variety of correction algorithms, such as Empirical Bayes (aka Combat), ANOVA and Median Polish.

Additional information can be found at http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview

# MBatch and MBatchUtils R Packages

The documentation directort contains several kinds of documentation for MBatch:

 * Files that start MBatch_01 are install documentations. Current instructions are for Linux (Debian 9.1). We expect to provide Windows and OS X instructions in late 2017/early 2018.
 * Files that start MBatch_02 are additional details about the test files in the package.
 * Files that start MBatch_03 are detail the file formats used by MBatch and the associated "Standardized Data" files.
 * Files that start MBatch_04 are documentation of assessment algorithms/plots.
 * Files that start MBatch_05 are documentation of correction algorithms.

Downloads and details on Standardized Data are available at http://bioinformatics.mdanderson.org/TCGA/databrowser/

## MBatch R Package

If you have the equivalent of Java 8 and R 3.4+ installed on your machine, and are familiar with your OS prerequisites and R package installation, the following quickstart instructions may allow quick installation.

```R
# required CRAN packages
install.packages(c("rJava", "devtools", "Cairo", "epiR", "gtools", "mclust", "squash", "httr"), dependencies=TRUE, repos = "http://cloud.r-project.org/")

# required GitHub package
library(devtools)
install_github("js229/Vennerable")

# required Bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite(c("limma","RBGL","graph","Biobase"), ask="a")
install.packages(c("oompaBase", "ClassDiscovery", "PreProcess"), dependencies=TRUE, repos=c("http://cloud.r-project.org", "http://silicovore.com/OOMPA/"))

## MBatch package
devtools::install_github("MD-Anderson-Bioinformatics/BatchEffectsPackage/apps/MBatch")
```

# MBatchUtils R Package

If you have the equivalent of Java 8 and R 3.4+ installed on your machine, and are familiar with your OS prerequisites and R package installation, the following quickstart instructions may allow quick installation.

First install MBatch, as provided above.

```R

# MBatchUtils package
library(devtools)
devtools::install_github("MD-Anderson-Bioinformatics/MBatchPackage/apps/MBatchUtils")
```


# MBatch Stand-Alone Docker Image with RStudio

The docker-compose.yml file provided with this release provides access to the Docker Hub image for MBatch, which includes RStudio. Connecting via 8080 using the username/password docker_rstudio/docker_rstudio allows access.

The docker-compose file uses the default host directories /MBatchSA/mbatch and /MBatchSA/website.
The docker-compose file can be started with: docker-compose -f docker-compose.yml start

Using "docker-compose -f docker-compose.yml start" will let you watch a "top" running in the guest Debian system.

Stop the container using: docker-compose -f docker-compose.yml down

Additional documentation in the yml file and coming soon.
