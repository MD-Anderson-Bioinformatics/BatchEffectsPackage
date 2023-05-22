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

message("starting installForMBatch1")

message("#*#* oompaBase, ClassDiscovery, and PreProcess")
# used by PCA+ and other clustering algorithms in MBatch
# uses ClassDiscovery functions distanceMatrix and SamplePCA
install.packages(c("oompaBase", "ClassDiscovery", "PreProcess"), dependencies=TRUE, repos=c("http://cran.r-project.org", "http://silicovore.com/OOMPA/"))

message("done installForMBatch1")
