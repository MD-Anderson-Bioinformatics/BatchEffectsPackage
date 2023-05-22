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

message("starting installForMBatch0")

message("#*#* reticulate used by MBatch")
install.packages("reticulate", dependencies=TRUE, repos = "http://cran.r-project.org")

message("#*#* mclust used by oompaBase, ClassDiscovery, and PreProcess (in install 1)")
install.packages("mclust", dependencies=TRUE, repos = "http://cran.r-project.org")

message("#*#* Biobase used by oompaBase, ClassDiscovery, and PreProcess (in install 1)")
install.packages("BiocManager", dependencies=TRUE, repos = "http://cran.r-project.org")
BiocManager::install(c("Biobase"), update=FALSE)

# no longer needed
#BiocManager::install(c("limma", "RBGL", "graph", "Biobase"), update=FALSE)

message("done installForMBatch0")

