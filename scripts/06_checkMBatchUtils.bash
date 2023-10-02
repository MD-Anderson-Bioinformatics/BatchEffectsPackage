#!/bin/bash

echo "START 06_checkMBatchUtils"

echo "06_checkMBatchUtils DSCCR_TEST_INPUT=${DSCCR_TEST_INPUT}"
echo "06_checkMBatchUtils DSCCR_TEST_OUTPUT=${DSCCR_TEST_OUTPUT}"

export DSCCR_TEST_INPUT=${DSCCR_TEST_INPUT}
export DSCCR_TEST_OUTPUT=${DSCCR_TEST_OUTPUT}

set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

export R_LIBS_USER="/tmp/R/4/lib:${R_LIBS_USER}"
# build/check does not use --library
env

# conda is a function, which is not propagated to bash scripts
# need to activate this so R can use it
echo "activate conda itself"
source /home/bcbuser/conda/etc/profile.d/conda.sh

echo "install PyMBatch"
cd PyMBatch
pip install .
cd ..

echo "build MBatch"
R CMD build MBatch

echo "install MBatch"
R CMD INSTALL MBatch_*.tar.gz

echo "build for check MBatchUtils"
R CMD build MBatchUtils

echo "list temp path"
ls -lh /tmp/R/4/lib

echo "install MBatchUtils"
R CMD INSTALL MBatchUtils_*.tar.gz

echo "build DscCR"
R CMD build DscCR
echo "check DscCR"
R CMD check DscCR*.tar.gz

echo "copy DscCR times"
SRC_FILE=/builds/BatchEffects_clean/BatchEffectsPackage/data/testing_dynamic/DSCCR/times.tsv
DST_DIR=/dsc_results
TIMESTAMP=$(date +"%Y_%m_%d_%H%M%S")
if [ -d "$DST_DIR" ];
then
    echo "$DST_DIR directory exists, copy output."
    if [ -f "$SRC_FILE" ];
    then
        echo "$SRC_FILE file exists, copy output."
        cp ${SRC_FILE} ${DST_DIR}/${TIMESTAMP}_times.tsv
    else
        echo "$SRC_FILE file does not exist, do not copy output."
        touch ${DST_DIR}/${TIMESTAMP}_times_NODATA.tsv
    fi
else
    echo "$DST_DIR directory does not exist."
fi

rm DscCR*.tar.gz

echo "FINISHED 06_checkMBatchUtils"

