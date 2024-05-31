# MBatch and MBatchUtils R Packages and MBatch Stand-Alone Docker Image

This is for educational and research purposes only. 

Samples from large research projects are often processed and run in multiple batches at different times. Because the samples are processed in batches rather than all at once, the data can be vulnerable to systematic noise such as batch effects (unwanted variation between batches) and trend effects (unwanted variation over time), which can lead to misleading analysis results.

The MBatch R package is designed to help assess and correct for batch effects. It first allows the user to assess and quantify the presence of any batch effects via algorithms such as Hierarchical Clustering, Principal Component Analysis, and box plots. If significant batch effects are observed in the data, the user then has the option of selecting from a variety of correction algorithms, such as Empirical Bayes (aka Combat), ANOVA and Median Polish.

Additional information can be found at http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview

Test data has been moved to a different repository, to address the issue where GitHub throttling code causes the devtools::install_github command to fail.
Test data can be added by running downloadData.bash, which is at the root of the repo.
The new data repo is https://github.com/MD-Anderson-Bioinformatics/BatchEffectsPackageData

|Component|Description|
|--|--|
|MBatch (R)|R package for basic Batch Effects|
|MBatchUtils (R)|R package for MutBatch, running from config, and making ZIP archives|
|BEVIndex (app)|Used by MBatchUtils to make ZIP Archives|

# MBatch Stand-Alone Docker Image

Download the docker-compose.yml file at the root of this repository. This file is setup for use on Linux.

Make the following directories.

 - /MBatchSA/mbatch
 - /MBatchSA/website

Permissions or ownership of the directories may need to be changed or matched to the Docker image user 2002.

In the directory with the docker-compose.yml file run:

	docker compose -f docker-compose.yml up --no-build -d

You can stop it with:

	docker compose -f docker-compose.yml down

To connect on the command line as the default user (2002) use:

	docker exec -t -i mbatchsa_cont_extr /bin/bash

See documentation at https://github.com/MD-Anderson-Bioinformatics/BatchEffectsPackage/blob/master/docs/MBatchImageSA/MBISA_01_InstallExtImageLinux.pdf for details.

# MBatch and MBatchUtils R Packages

# **Install documentation is in the process of being updated.**

The documentation directort contains several kinds of documentation for MBatch:

 * Files that start MBatch_01 are install documentations.
 * Files that start MBatch_02 are additional details about the test files in the package.
 * Files that start MBatch_03 are detail the file formats used by MBatch and the associated "Standardized Data" files.
 * Files that start MBatch_04 are documentation of assessment algorithms/plots.
 * Files that start MBatch_05 are documentation of correction algorithms.

Downloads and details on Standardized Data are available at http://bioinformatics.mdanderson.org/TCGA/databrowser/

## MBatch Python Environment

MBatch requires a Python environment with Anaconda.

Something similar to this should allow setup on Linux.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
mkdir /home/bcbuser/conda
# do not use unattended install so you can select automatic init
bash /home/bcbuser/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh -p /home/bcbuser/conda -f 
source /home/bcbuser/conda/bin/activate
conda init
conda update -y conda
```

Then do required installs, similar to this:

```
conda create -y -n gendev
conda activate gendev
conda install -y -c conda-forge python==3.11.*
conda install -y -c conda-forge pandas
conda install -y -c conda-forge numpy
conda install -y -c conda-forge matplotlib
conda install -y -c conda-forge scanpy
conda install -y -c conda-forge pillow
conda install -y -c conda-forge jsonpickle
conda install -y -c conda-forge requests
conda install -y -c conda-forge xmltodict
conda install -y -c conda-forge cryptography
conda install -y -c conda-forge urllib3
conda install -y -c conda-forge scipy
```
Then you can install the MBatch Python Package:

```
conda activate gendev
pip install "git+https://github.com/MD-Anderson-Bioinformatics/BatchEffectsPackage.git#egg=mbatch&subdirectory=apps/PyMBatch"
```

If needed, you must override the default Python environment of "/BEA/gendev" by setting the environmental variable MBATCH_PYTHON_ENV prior to install. On Linux or OS X, you can create a link from your gendev to /BEA/gendev.

## MBatch R Package

If you are familiar with your OS prerequisites and R package installation, the following quickstart instructions may allow quick installation.

```R
install.packages(c("BiasedUrn", "Cairo", "covr", "devtools", "dunn.test", "gert", "htmlwidgets", "httr", "jsonlite", "lubridate", "magick", "mclust", "pander", "reticulate", "rversions", "sf", "shiny", "squash", "usethis", "uwot"), dependencies=TRUE, repos = "http://cran.r-project.org")

#Install version 2.0.40, anything newer fails to install
install.packages("https://cran.r-project.org/src/contrib/Archive/epiR/epiR_2.0.40.tar.gz", dependencies=TRUE, repos = "http://cran.r-project.org")

install.packages("BiocManager", dependencies=TRUE, repos = "http://cran.r-project.org")
BiocManager::install(c("Biobase"), update=FALSE)

install.packages(c("oompaBase", "ClassDiscovery", "PreProcess"), dependencies=TRUE, repos=c("http://cran.r-project.org", "http://silicovore.com/OOMPA/"))

library(devtools)
devtools::install_github("jlmelville/vizier")
devtools::install_github('MD-Anderson-Bioinformatics/tsvio')
devtools::install_github('MD-Anderson-Bioinformatics/NGCHMSupportFiles', ref='main')
devtools::install_github('MD-Anderson-Bioinformatics/NGCHM-R')

## MBatch package
devtools::install_github("MD-Anderson-Bioinformatics/BatchEffectsPackage/apps/MBatch")
```
Test data has been moved to a different repository, to address the issue where GitHub throttling code causes the devtools::install_github command to fail. See the Linux install documentation on using test data.

**For educational and research purposes only.**

**Funding** 
This work was supported in part by U.S. National Cancer Institute (NCI) grant: Weinstein, Mills, Akbani. Batch effects in molecular profiling data on cancers: detection, quantification, interpretation, and correction, 5U24CA210949

