#!/bin/bash

echo "START 03_checkMBatch"
set -e
BASE_DIR=$1

echo "move to app dir"
cd ${BASE_DIR}/apps

# conda is a function, which is not propagated to bash scripts
# need to activate this so R can use it
echo "activate conda itself"
source /home/bcbuser/conda/etc/profile.d/conda.sh

echo "cd PyMBatch"
cd PyMBatch
echo "conda activate gendev"
conda activate gendev
echo "install PyMBatch for MBatch R build"
pip install .

echo "build for check MBatch"
cd ..
R CMD build MBatch

echo "check MBatch"
env _R_CHECK_FORCE_SUGGESTS_=0 R CMD check MBatch_*.tar.gz

echo "remove build"
rm MBatch_*.tar.gz

echo "FINISHED 03_checkMBatch"

