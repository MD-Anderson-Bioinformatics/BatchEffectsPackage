#!/bin/bash

echo "START 08_copyForDocker"
set -e
BASE_DIR=$1

echo "remove and create installations directory"
# rm with Force, to not fail if it doesn't exist (which it won't after a git clone)
# then mkdir
rm -rf ${BASE_DIR}/docker-build/MBatchImage/installations
mkdir ${BASE_DIR}/docker-build/MBatchImage/installations

# do not do this here -- is too big for default artifact storage
#echo "copy testing_static"
#cp -r ${BASE_DIR}/data/testing_static ${BASE_DIR}/docker-build/MBatchImage/installations/.

echo "copy R packages"
cp ${BASE_DIR}/apps/*.tar.gz ${BASE_DIR}/docker-build/MBatchImage/installations/.

echo "list files"
ls -lh ${BASE_DIR}/docker-build/MBatchImage/installations/*

echo "FINISHED 08_copyForDocker"

