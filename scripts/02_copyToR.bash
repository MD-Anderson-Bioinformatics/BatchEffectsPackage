#!/bin/bash

echo "START 02_copyToR"
set -e
BASE_DIR=$1

echo "copy DscJava to MBatch"
ls -lh ${BASE_DIR}/apps/DscJava/target
rm -rf ${BASE_DIR}/apps/MBatch/inst/DscJava
ls -lh ${BASE_DIR}/apps/MBatch/inst
mkdir ${BASE_DIR}/apps/MBatch/inst/DscJava
cp ${BASE_DIR}/apps/DscJava/target/*.jar ${BASE_DIR}/apps/MBatch/inst/DscJava/.
mv ${BASE_DIR}/apps/MBatch/inst/DscJava/DscJava-*.jar ${BASE_DIR}/apps/MBatch/inst/DscJava/DscJava.jar
ls -lh ${BASE_DIR}/apps/MBatch/inst/DscJava

echo "copy LegendJava to MBatch"
ls -lh ${BASE_DIR}/apps/LegendJava/target
rm -rf ${BASE_DIR}/apps/MBatch/inst/LegendJava
mkdir ${BASE_DIR}/apps/MBatch/inst/LegendJava
cp ${BASE_DIR}/apps/LegendJava/target/*.jar ${BASE_DIR}/apps/MBatch/inst/LegendJava/.
mv ${BASE_DIR}/apps/MBatch/inst/LegendJava/LegendJava-*.jar ${BASE_DIR}/apps/MBatch/inst/LegendJava/LegendJava.jar
ls -lh ${BASE_DIR}/apps/MBatch/inst/LegendJava

echo "copy ReadRJava to MBatch"
ls -lh ${BASE_DIR}/apps/ReadRJava/target
rm -rf ${BASE_DIR}/apps/MBatch/inst/ReadRJava
mkdir ${BASE_DIR}/apps/MBatch/inst/ReadRJava
ls -lh ${BASE_DIR}/apps/MBatch/inst/ReadRJava
cp ${BASE_DIR}/apps/ReadRJava/target/*.jar ${BASE_DIR}/apps/MBatch/inst/ReadRJava/.
mv ${BASE_DIR}/apps/MBatch/inst/ReadRJava/ReadRJava-*.jar ${BASE_DIR}/apps/MBatch/inst/ReadRJava/ReadRJava.jar
ls -lh ${BASE_DIR}/apps/MBatch/inst/ReadRJava

echo "copy BEVIndex to MBatchUtils"
ls -lh ${BASE_DIR}/apps/BEVIndex/target
rm -rf ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex
mkdir ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex
ls -lh ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex
cp ${BASE_DIR}/apps/BEVIndex/target/*.jar ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex/.
mv ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex/BEVIndex-*.jar ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex/BEVIndex.jar
ls -lh ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex

echo "copy tooltip files for BEVIndex to MBatchUtils"
cp ${BASE_DIR}/data/testing_static/BEVIndex/data/* ${BASE_DIR}/apps/MBatchUtils/inst/BEVIndex/.

echo "list files"
ls -lh ${BASE_DIR}/apps/MBatch/inst/*
ls -lh ${BASE_DIR}/apps/MBatchUtils/inst/*

echo "FINISHED 02_copyToR"

