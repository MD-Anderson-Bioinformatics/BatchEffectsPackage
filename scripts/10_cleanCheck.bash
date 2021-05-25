#!/bin/bash

# not done during CI
echo "START 10_cleanCheck"
set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

echo "remove .Rcheck directories"
rm -rf *.Rcheck

echo "remove R packages"
rm -f *.tar.gz

echo "remove Docker installations data"
rm -rf ${BASE_DIR}/docker-build/MBatchImage/installations/*

echo "remove JAR files in R packages"
rm -rf ${BASE_DIR}/apps/MBatch/inst/DscJava/*
rm -rf ${BASE_DIR}/apps/MBatch/inst/LegendJava/*
rm -rf ${BASE_DIR}/apps/MBatch/inst/ReadRJava/*
rm -rf ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex/*

echo "list apps directory"
ls -lh ${BASE_DIR}/apps

echo "list installations directory"
ls -lh ${BASE_DIR}/docker-build/MBatchImage/installations

echo "FINISHED 10_cleanCheck"

