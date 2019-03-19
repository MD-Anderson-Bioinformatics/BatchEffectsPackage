#!/bin/bash

echo "Clear build and dist files"
rm -r ../apps/*/build/*
rm -r ../apps/*/dist/*

echo "Clear test output data"
rm -r ../data/testing_dynamic/BEVIndex/*
rm -r ../data/testing_dynamic/Boxplot/*
rm -r ../data/testing_dynamic/LegendJava/*
rm -r ../data/testing_dynamic/MBatch/*
rm -r ../data/testing_dynamic/MBatchUtils/*
rm -r ../data/testing_dynamic/ReadRJava/*

echo "ls build and dist files"
ls ../apps/*/build
ls ../apps/*/dist

echo "ls old output data"
ls -l ../data/testing_dynamic/BEVIndex
ls -l ../data/testing_dynamic/Boxplot
ls -l ../data/testing_dynamic/LegendJava
ls -l ../data/testing_dynamic/MBatch
ls -l ../data/testing_dynamic/MBatchUtils
ls -l ../data/testing_dynamic/ReadRJava

echo "done"
