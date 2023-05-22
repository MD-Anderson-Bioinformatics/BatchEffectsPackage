# MBatch Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>


GLOBAL_MBATCH_ENV = new.env(parent = emptyenv())

####
#### variable to hold optional TEST version (for directory name)
####

assign('MBATCH_ERROR_TEST', NULL, GLOBAL_MBATCH_ENV)

getGlobalMBatchErrorTest <- function()
{
  return(get('MBATCH_ERROR_TEST', envir = GLOBAL_MBATCH_ENV))
}

setGlobalMBatchErrorTest <- function(theValue)
{
  assign('MBATCH_ERROR_TEST', theValue, GLOBAL_MBATCH_ENV)
  return(getGlobalMBatchErrorTest())
}

####
#### variable to hold optional TEST version (for directory name)
####

assign('MBATCH_VERSION_TEST', NULL, GLOBAL_MBATCH_ENV)

getGlobalMBatchVersionTest <- function()
{
  return(get('MBATCH_VERSION_TEST', envir = GLOBAL_MBATCH_ENV))
}

setGlobalMBatchVersionTest <- function(theValue)
{
  assign('MBATCH_VERSION_TEST', theValue, GLOBAL_MBATCH_ENV)
  return(getGlobalMBatchVersionTest())
}

####
#### variable to hold optional DATA version (for directory name)
####

assign('MBATCH_VERSION_DATA', NULL, GLOBAL_MBATCH_ENV)

getGlobalMBatchVersionData <- function()
{
  return(get('MBATCH_VERSION_DATA', envir = GLOBAL_MBATCH_ENV))
}

setGlobalMBatchVersionData <- function(theValue)
{
  assign('MBATCH_VERSION_DATA', theValue, GLOBAL_MBATCH_ENV)
  return(getGlobalMBatchVersionData())
}

####
#### variable to hold Python environment name
####

assign('MBATCH_PYTHON_ENV', "gendev", GLOBAL_MBATCH_ENV)

getGlobalMBatchEnv <- function()
{
  return(get('MBATCH_PYTHON_ENV', envir = GLOBAL_MBATCH_ENV))
}

setGlobalMBatchEnv <- function(theValue)
{
  assign('MBATCH_PYTHON_ENV', theValue, GLOBAL_MBATCH_ENV)
  use_condaenv(getGlobalMBatchEnv())
  return(getGlobalMBatchEnv())
}

####
#### variable to hold logging information
####

assign('MBATCH_LOGGER', new("Logging"), GLOBAL_MBATCH_ENV)

getGlobalMBatchLogger <- function()
{
  return(get('MBATCH_LOGGER', envir = GLOBAL_MBATCH_ENV))
}

setGlobalMBatchLogger <- function(theValue)
{
  assign('MBATCH_LOGGER', theValue, GLOBAL_MBATCH_ENV)
  return(getGlobalMBatchLogger())
}

