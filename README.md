# MBatch and MBatchUtils R Packages and MBatch Stand-Alone Docker Image

This is for educational and research purposes only. 

Samples from large research projects are often processed and run in multiple batches at different times. Because the samples are processed in batches rather than all at once, the data can be vulnerable to systematic noise such as batch effects (unwanted variation between batches) and trend effects (unwanted variation over time), which can lead to misleading analysis results.

The MBatch R package is designed to help assess and correct for batch effects. It first allows the user to assess and quantify the presence of any batch effects via algorithms such as Hierarchical Clustering, Principal Component Analysis, and box plots. If significant batch effects are observed in the data, the user then has the option of selecting from a variety of correction algorithms, such as Empirical Bayes (aka Combat), ANOVA and Median Polish.

Additional information can be found at http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview

|Component|Description|
|--|--|
|MBatch (R)|R package for basic Batch Effects|
|MBatchUtils (R)|R package for MutBatch, running from config, and making ZIP archives|
|BEVIndex (app)|Used by MBatchUtils to make ZIP Archives|
|DscJava (app)|Java code for computing DSC values, used by MBatch|
|LegendJava (app)|Java code for writing legends, used by MBatch|
|ReadRJava (app)|Java code for reading large TSVs, used by MBatch|
|DebianRJava (docker)|Base image with R and Java|
|MBatchImageSA (docker)|Image with MBatch R packages|

# MBatch Stand-Alone Docker Image

Download the docker-compose.yml file at the root of this repository. This file is setup for use on Linux.

Make the following directories.

 - /MBatchSA/mbatch
 - /MBatchSA/website

Permissions or ownership of the directories may need to be changed or matched to the Docker image user 2002.

In the directory with the docker-compose.yml file run:

	docker-compose -f docker-compose.yml up --no-build -d

You can stop it with:

	docker-compose -f docker-compose.yml down

To connect on the command line as the default user (2002) use:

	docker exec -t -i mbatchsa_cont_extr /bin/bash

See documentation at https://github.com/MD-Anderson-Bioinformatics/BatchEffectsPackage/blob/master/docs/MBatchImageSA/MBISA_01_InstallExtImageLinux.pdf for details.

# MBatch and MBatchUtils R Packages

The documentation directort contains several kinds of documentation for MBatch:

 * Files that start MBatch_01 are install documentations.
 * Files that start MBatch_02 are additional details about the test files in the package.
 * Files that start MBatch_03 are detail the file formats used by MBatch and the associated "Standardized Data" files.
 * Files that start MBatch_04 are documentation of assessment algorithms/plots.
 * Files that start MBatch_05 are documentation of correction algorithms.

Downloads and details on Standardized Data are available at http://bioinformatics.mdanderson.org/TCGA/databrowser/

## MBatch R Package

If you have the equivalent of Java 8 and R 4+ installed on your machine, and are familiar with your OS prerequisites and R package installation, the following quickstart instructions may allow quick installation.

```R
# required CRAN packages
install.packages(c("rJava", "devtools", "Cairo", "epiR", "gtools", "mclust", "squash", "httr", "eulerr"), dependencies=TRUE, repos = "http://cloud.r-project.org/")

# required Bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite(c("limma","RBGL","graph","Biobase"), ask="a")
install.packages(c("oompaBase", "ClassDiscovery", "PreProcess"), dependencies=TRUE, repos=c("http://cloud.r-project.org", "http://silicovore.com/OOMPA/"))

## MBatch package
devtools::install_github("MD-Anderson-Bioinformatics/BatchEffectsPackage/apps/MBatch")
```

# MBatchUtils R Package

If you have the equivalent of Java 8 and R 4+ installed on your machine, and are familiar with your OS prerequisites and R package installation, the following quickstart instructions may allow quick installation.

First install MBatch, as provided above.

```R

# MBatchUtils package
library(devtools)
devtools::install_github("MD-Anderson-Bioinformatics/MBatchPackage/apps/MBatchUtils")
```

**For educational and research purposes only.**

**Funding** 
This work was supported in part by U.S. National Cancer Institute (NCI) grant: Weinstein, Mills, Akbani. Batch effects in molecular profiling data on cancers: detection, quantification, interpretation, and correction, 5U24CA210949

