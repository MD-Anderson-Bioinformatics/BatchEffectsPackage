#!/bin/bash

echo "START 02_copyToR"
set -e
BASE_DIR=$1

echo "Link/populate /BEA/BatchEffectsPackage_data/testing_static"

echo "ln -s ${BASE_DIR}/data/testing_static /BEA/BatchEffectsPackage_data/testing_static"

mkdir /BEA/BatchEffectsPackage_data
ln -s ${BASE_DIR}/data/testing_static /BEA/BatchEffectsPackage_data/testing_static

echo "FINISHED 02_copyToR"

