FROM mdabcb/bcb_base:2024-01-11-0900

# NOTE: above version should match build from Dockerfile for BCBBase, which is manually done.
#       manual build of BCBBase is because of the extended build time to generate the image.
#       it is generally rebuilt very seldomly.

## This can be built with: docker build -t mdabcb/mbatchsa_image:HUB_2024-05-29-1000 .
## Login to docker hub with: docker login
## repo may need be created on Docker Hub before pushing (depends on versions of stuff)
## Push to Docker Hub with: docker push mdabcb/mbatchsa_image:HUB_2024-05-29-1000

# reminder, use {} around environmental variables, otherwise docker uses it as a literal

LABEL edu.mda.bcb.name="MBatch" \
      edu.mda.bcb.sub="mba" \
      edu.mda.bcb.mba.version="2024-05-29-1000"

# make bcbuser
RUN groupadd -g 2003 bcbgroup && \
    useradd -l -s /bin/bash -d /home/bcbuser -u 2002 -g bcb_base_group bcbuser && \
    usermod -a -G bcbgroup bcbuser && \
    chown bcbuser:bcbgroup /home/bcbuser && \
    mkdir /home/bcbuser/mbatch && \
    chown -R :bcb_base_group /home/bcbuser/mbatch && \
    mkdir /BEA/BatchEffectsPackage_data && \
    mkdir /BEA/BatchEffectsPackage_data/testing_static && \
    chown -R :bcb_base_group /BEA/BatchEffectsPackage_data && \
    chmod -R 775 /BEA/BatchEffectsPackage_data

# COPY is always done as root!!!!
# quote: When using COPY with more than one source file, the destination must be a directory and end with a / :unquote
COPY --chown=bcbuser:bcbgroup --chmod=775 installations/*.tar.gz /home/bcbuser/mbatch/
#################################################
# this COPY block adds 2.8GB of layer data
#################################################
COPY --chown=bcbuser:bcbgroup --chmod=775 installations/testing_static /BEA/BatchEffectsPackage_data/testing_static
COPY --chown=bcbuser:bcbgroup --chmod=775 installations/PyMBatch /home/bcbuser/PyMBatch/
# perform Python packge install
RUN . /home/bcbuser/conda/bin/activate && \
    conda activate /home/bcbuser/conda/envs/gendev && \
    cd /home/bcbuser/PyMBatch && \
    pip install .

# perform R installs for packages, and install MBatch and MBatchUtils
RUN cd /home/bcbuser/mbatch && \
    R CMD INSTALL MBatch_*.tar.gz && \
    R CMD INSTALL MBatchUtils_*.tar.gz

ENV PATH="/home/bcbuser/conda/bin:$PATH"

