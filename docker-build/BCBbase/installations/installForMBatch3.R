# Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

message("starting installForMBatch3")

# used in DAPIR
message("#*#* httr")
install.packages("httr", dependencies=TRUE, repos = "http://cran.r-project.org")
message("#*#* jsonlite")
install.packages("jsonlite", dependencies=TRUE, repos = "http://cran.r-project.org")

# uses dunn.test in MBatchUtils
message("#*#* dunn.test")
install.packages("dunn.test", dependencies=TRUE, repos = "http://cran.r-project.org")

# need timeout for larger packages
options(timeout=9999999)

# used for UMAP
message("#*#* uwot")
install.packages("uwot", dependencies=TRUE, repos = "http://cran.r-project.org")

# required for devtools (UMAP install and elsewhere)
message("#*#* stringi")
install.packages("stringi", dependencies=TRUE, repos = "http://cran.r-project.org")
message("#*#* devtools")
install.packages("devtools", dependencies=TRUE, repos = "http://cran.r-project.org")
library(devtools)

# required for UMAP embeds
message("#*#* jlmelville/vizier")
devtools::install_github("jlmelville/vizier")

message("done installForMBatch3")
