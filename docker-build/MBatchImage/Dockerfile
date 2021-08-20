FROM mdabcb/bcb_base:2021-08-13-1100

## This can be built with: docker build -t mdabcb/mbatchsa_image:BEA_VERSION_TIMESTAMP .
## Login to docker hub with: docker login
## repo may need be created on Docker Hub before pushing (depends on versions of stuff)
## Push to Docker Hub with: docker push mdabcb/mbatchsa_image:BEA_VERSION_TIMESTAMP

# reminder, use {} around environmental variables, otherwise docker uses it as a literal

LABEL edu.mda.bcb.name="MBatch" \
      edu.mda.bcb.sub="mba" \
      edu.mda.bcb.mba.version="BEA_VERSION_TIMESTAMP"

# make bcbuser
RUN mkdir /home/bcbuser && \
    mkdir /home/bcbuser/mbatch && \
    useradd -l -s /bin/bash -d /home/bcbuser -u 2002 bcbuser && \
    echo "bcbuser:bcbuser" | chpasswd && \
    chown -R bcbuser:bcbuser /home/bcbuser && \
    mkdir /BatchEffectsPackage_data && \
    mkdir /BatchEffectsPackage_data/testing_static && \
    chown -R bcbuser:bcbuser /BatchEffectsPackage_data && \
    chmod -R u+rwx /BatchEffectsPackage_data

# COPY is always done as root!!!!
# quote: When using COPY with more than one source file, the destination must be a directory and end with a / :unquote
COPY installations/*.tar.gz /home/bcbuser/mbatch/
COPY installations/testing_static /BatchEffectsPackage_data/testing_static
RUN ls -lh /BatchEffectsPackage_data && \
    chown -R bcbuser:bcbuser /BatchEffectsPackage_data && \
    chmod -R u+rwx /BatchEffectsPackage_data && \
    mkdir /bea_testing && \
    ln -s /BatchEffectsPackage_data/testing_static/MATRIX_DATA /bea_testing/MATRIX_DATA && \
    chown -R bcbuser:bcbuser /bea_testing && \
    chmod -R u+rwx /bea_testing

# perform R installs for packages, and install MBatch and MBatchUtils
RUN cd /home/bcbuser/mbatch && \
    R CMD INSTALL MBatch_*.tar.gz && \
    R CMD INSTALL MBatchUtils_*.tar.gz

