#!/bin/bash

echo "START 03_checkMBatch"
set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

echo "build for check MBatch"
R CMD build MBatch

echo "check MBatch"
env _R_CHECK_FORCE_SUGGESTS_=0 R CMD check MBatch_*.tar.gz

echo "remove build"
rm MBatch_*.tar.gz

echo "FINISHED 03_checkMBatch"

