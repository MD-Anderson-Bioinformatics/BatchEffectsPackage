#!/bin/bash

echo "START 07_buildMBatchUtils"
set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

echo "list temp path"
ls -lh /tmp/R/4/lib

# conda is a function, which is not propagated to bash scripts
# need to activate this so R can use it
echo "activate conda itself"
source /home/bcbuser/conda/etc/profile.d/conda.sh

echo "build MBatchUtils"
export R_LIBS_USER="/tmp/R/4/lib:${R_LIBS_USER}"
R CMD build MBatchUtils

echo "FINISHED 07_buildMBatchUtils"

