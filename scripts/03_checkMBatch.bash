#!/bin/bash

echo "START 03_checkMBatch"
set -e
BASE_DIR=$1

# conda is a function, which is not propagated to bash scripts
# need to activate this so R can use it
echo "activate conda itself"
source /home/bcbuser/conda/etc/profile.d/conda.sh
echo "conda activate gendev"
conda activate gendev
echo "python version"
python --version

echo "move to app dir"
cd ${BASE_DIR}/apps

echo "cd PyMBatch"
cd PyMBatch
echo "install PyMBatch for MBatch R build"
pip install .

echo "build for check MBatch"
cd ..
R CMD build MBatch

echo "check MBatch"
env _R_CHECK_FORCE_SUGGESTS_=0 R CMD check MBatch_*.tar.gz

echo "******************************"
echo "******************************"
cat ${BASE_DIR}/apps/MBatch.Rcheck/00check.log
echo "******************************"
echo "******************************"

echo "remove build"
rm MBatch_*.tar.gz

echo "FINISHED 03_checkMBatch"

