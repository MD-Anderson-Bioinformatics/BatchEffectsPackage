#!/bin/bash

echo "START 04_buildMBatch"
set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

# conda is a function, which is not propagated to bash scripts
# need to activate this so R can use it
echo "activate conda itself"
source /home/bcbuser/conda/etc/profile.d/conda.sh

echo "build MBatch"
R CMD build MBatch

echo "make temp install dir"
mkdir -p /tmp/R/4/lib

echo "install MBatch -- for use by MBatchUtils"
R CMD INSTALL --library=/tmp/R/4/lib MBatch_*.tar.gz

echo "FINISHED 04_buildMBatch"

