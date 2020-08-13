#!/bin/bash

echo "move to app dir"
cd ../apps

echo "prepare copy to MBatch"
cd MBatch
cd inst

echo "copy DscJava to MBatch"
cp ../../DscJava/dist/*.jar DscJava/.
cp ../../DscJava/dist/lib/*.jar DscJava/.

echo "copy LegendJava to MBatch"
cp ../../LegendJava/dist/*.jar LegendJava/.
cp ../../LegendJava/dist/lib/*.jar LegendJava/.

echo "copy ReadRJava to MBatch"
cp ../../ReadRJava/dist/*.jar ReadRJava/.
cp ../../ReadRJava/dist/lib/*.jar ReadRJava/.

echo "prepare copy to MBatchUtils"
cd ../..
cd MBatchUtils
cd inst

echo "copy BEVIndex to MBatchUtils"
cp ../../BEVIndex/dist/*.jar BEVIndex/.
cp ../../BEVIndex/dist/lib/*.jar BEVIndex/.

echo "copy tooltip files for BEVIndex to MBatchUtils"
cp ../../../data/testing_static/BEVIndex/data/* BEVIndex/.

echo "list files"
cd ../..
ls -l MBatch/inst/*
ls -l MBatchUtils/inst/*

echo "done"
