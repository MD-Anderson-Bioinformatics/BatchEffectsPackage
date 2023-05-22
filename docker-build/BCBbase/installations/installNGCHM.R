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

message("starting")

message("#*#* httr")
install.packages("httr", dependencies = TRUE, repos = "http://cran.r-project.org")
message("#*#* magick")
install.packages("magick", dependencies = TRUE, repos = "http://cran.r-project.org")
message("#*#* gert")
install.packages("gert", dependencies = TRUE, repos = "http://cran.r-project.org")
message("#*#* htmlwidgets - do not set dependecies. Is circular relationship with shiny")
install.packages("htmlwidgets", dependencies = FALSE, repos = "http://cran.r-project.org")
message("#*#* shiny")
install.packages("shiny", dependencies = TRUE, repos = "http://cran.r-project.org")
message("#*#* usethis")
install.packages("usethis", dependencies = TRUE, repos = "http://cran.r-project.org")
message("#*#* covr")
install.packages("covr", dependencies = TRUE, repos = "http://cran.r-project.org")
message("#*#* rversions")
install.packages("rversions", dependencies = TRUE, repos = "http://cran.r-project.org")
message("#*#* devtools")
install.packages("devtools", dependencies = TRUE, repos = "http://cran.r-project.org")
library(devtools)
# remotes
message("#*#* MD-Anderson-Bioinformatics/tsvio")
devtools::install_github('MD-Anderson-Bioinformatics/tsvio')
message("#*#* MD-Anderson-Bioinformatics/NGCHMSupportFiles")
devtools::install_github('MD-Anderson-Bioinformatics/NGCHMSupportFiles', ref='main')
message("#*#* MD-Anderson-Bioinformatics/NGCHM-R")
devtools::install_github('MD-Anderson-Bioinformatics/NGCHM-R')

message("done")
