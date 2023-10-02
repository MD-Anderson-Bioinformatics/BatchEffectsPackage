# PyMBatch

## MBatch Python Environment

MBatch requires a Python environment with Anaconda.

Something similar to this should allow setup on Linux.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
mkdir /home/bcbuser/conda
# do not use unattended install so you can select automatic init
bash /home/bcbuser/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh -p /home/bcbuser/conda -f 
source /home/bcbuser/conda/bin/activate
conda init
conda update -y conda
```

Then do required installs, similar to this:

```
conda create -y -n gendev
conda activate gendev
conda install -y -c conda-forge python==3.10
conda install -y -c conda-forge pandas
conda install -y -c conda-forge numpy
conda install -y -c conda-forge matplotlib
conda install -y -c conda-forge pillow
conda install -y -c conda-forge jsonpickle
conda install -y -c conda-forge xmltodict
conda install -y -c conda-forge pillow
```
Then you can install the MBatch Python Package:

```
conda activate gendev
pip install git+https://github.com/MD-Anderson-Bioinformatics/BatchEffectsPackage.git#egg=mbatch&subdirectory=apps/PyMBatch
```

