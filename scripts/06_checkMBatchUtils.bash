#!/bin/bash

echo "START 06_checkMBatchUtils"
set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

export R_LIBS_USER="/tmp/R/4/lib:${R_LIBS_USER}"
# build/check does not use --library
env

echo "build for check MBatchUtils"
R CMD build MBatchUtils

echo "list temp path"
ls -lh /tmp/R/4/lib

echo "check MBatchUtils"
R CMD check MBatchUtils_*.tar.gz

echo "remove build"
rm MBatchUtils_*.tar.gz

echo "FINISHED 06_checkMBatchUtils"

