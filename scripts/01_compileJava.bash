#!/bin/bash

echo "START 01_compileJava"
set -e
BASE_DIR=$1

echo "compile BEVIndex"
cd ${BASE_DIR}/apps/BEVIndex
mvn clean install dependency:copy-dependencies

echo "compile DscJava"
cd ${BASE_DIR}/apps/DscJava
mvn clean install dependency:copy-dependencies

echo "compile LegendJava"
cd ${BASE_DIR}/apps/LegendJava
mvn clean install dependency:copy-dependencies

echo "compile ReadRJava"
cd ${BASE_DIR}/apps/ReadRJava
mvn clean install dependency:copy-dependencies

echo "list target files"
ls -l ${BASE_DIR}/apps/*/target/*.jar

# turn off error checking, since rm/find generates a harmless error
set +e
echo "clean maven dirs"
find ${BASE_DIR} -name "\?" -type d -exec rm -rf {} \;

echo "FINISHED 01_compileJava"

