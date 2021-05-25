#!/bin/bash

echo "START 07_buildMBatchUtils"
set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

echo "list temp path"
ls -lh /tmp/R/4/lib

echo "build MBatchUtils"
export R_LIBS_USER="/tmp/R/4/lib:${R_LIBS_USER}"
R CMD build MBatchUtils

echo "FINISHED 07_buildMBatchUtils"

