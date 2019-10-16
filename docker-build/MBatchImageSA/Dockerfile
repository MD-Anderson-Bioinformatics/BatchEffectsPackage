FROM mdabcb/debian_r_java:2019-06-24-1600

# reminder, use {} around environmental variables, otherwise docker uses it as a literal

LABEL edu.mda.bcb.name="MBatch Stand Alone" \
      edu.mda.bcb.version="2019-10-15-1100"

# create and setup docker_rstudio
# also, set permissions and ownerships on internal docker directories
RUN mkdir /home/docker_rstudio && \
    useradd -l -s /bin/bash -d /home/docker_rstudio -u 2002 docker_rstudio && \
    echo "docker_rstudio:docker_rstudio" | chpasswd && \
    chown -R docker_rstudio:docker_rstudio /home/docker_rstudio

# create and setup docker_rstudio user
# also, set permissions and ownerships on internal docker directories
# make and set ownership and permissions on internal image volume mount locations
# must be done in one step, or could lock permissions in future images
RUN mkdir -p /MBATCH && \
    chown -R docker_rstudio:docker_rstudio /MBATCH && \
    chmod -R u+rwx /MBATCH && \
    mkdir -p /WEBSITE && \
    chown -R docker_rstudio:docker_rstudio /WEBSITE && \
    chmod -R u+rwx /WEBSITE

# make subdirectory in user home for running
# COPY is always done as root!!!!
RUN mkdir /home/docker_rstudio/mbatch && \
    chown -R docker_rstudio:docker_rstudio /home/docker_rstudio/mbatch && \
    chmod -R u+rwx /home/docker_rstudio/mbatch
COPY installations /home/docker_rstudio/mbatch

#update for RStudio Server install packages
RUN apt-get update && \
    apt-get install -f -y libapparmor1 libclang-3.8-dev libclang-common-3.8-dev libclang-dev libclang1-3.8 libgc1c2 libllvm3.8 libobjc-6-dev libobjc4 sudo && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install R studio
RUN chown -R docker_rstudio:docker_rstudio /home/docker_rstudio && \
    chmod -R u+rwx /home/docker_rstudio && \
    cd /home/docker_rstudio/mbatch && \
    ls -l && \
    gdebi -n rstudio-server-1.2.5001-amd64.deb

# perform R installs for packages, and install MBatch and MBatchUtils
RUN cd /home/docker_rstudio/mbatch && \
    R CMD INSTALL --install-tests MBatch_*.tar.gz && \
    R CMD INSTALL --install-tests MBatchUtils_*.tar.gz && \
    rm /home/docker_rstudio/mbatch/*.tar.gz && \
    rm /home/docker_rstudio/mbatch/*.deb

# set working directory for future commands
WORKDIR /home/docker_rstudio/mbatch
CMD ["/bin/bash", "./runRstudioServer.bash"]
