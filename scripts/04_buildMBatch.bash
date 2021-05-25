#!/bin/bash

echo "START 04_buildMBatch"
set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

echo "build MBatch"
R CMD build MBatch

echo "make temp install dir"
mkdir -p /tmp/R/4/lib

echo "install MBatch -- for use by MBatchUtils"
R CMD INSTALL --library=/tmp/R/4/lib MBatch_*.tar.gz

echo "FINISHED 04_buildMBatch"

