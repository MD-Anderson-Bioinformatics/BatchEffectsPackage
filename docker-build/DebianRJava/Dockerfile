FROM debian:9

# reminder, use {} around environmental variables, otherwise docker uses it as a literal

LABEL edu.mda.bcb.name="DebianRJava" \
      edu.mda.bcb.version="2020-09-11-1000" \
      edu.mda.bcb.Rversion="3.6" \
      edu.mda.bcb.Javaversion="8" \
      edu.mda.bcb.Debianversion="9"

# install JDK, apt utils
# separate layers to create more/smaller layers, to attempt to fix pull issue
RUN apt-get update && \
    apt-get upgrade -f -y && \
    apt-get install -f -y apt-utils readline-common procps curl && \
    apt-get install -f -y --no-install-recommends openjdk-8-jdk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# not sure why, but this does not work as part of the other apt-get chain
RUN apt-get update && \
    apt-get install -f -y gnupg && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# R requirements
RUN apt-get update && \
    apt-get install -f -y texlive && \
    apt-get install -f -y wget libssh2-1 libcairo2-dev libxt-dev curl && \
    apt-get install -f -y wget libcurl4-openssl-dev libssl1.0-dev libxml2-dev && \
    apt-get install -f -y git openssh-client && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Installing R 3.5+ on Debian https://cran.r-project.org/bin/linux/debian/
# do apt-get update and install, if this fails with "temporary error", firewall may be blocking DNS, check your docker DNS entries
# https://stackoverflow.com/questions/24991136/docker-build-could-not-resolve-archive-ubuntu-com-apt-get-fails-to-install-a?rq=1
# apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
# apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
####    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
# The above is constantly failing for one reason or another. Use saved key in file
####
# gdebi-core is used for installing RStudio later
COPY installations/jranke.key /jranke.key

# use newer version of R, currently 3.6
RUN echo 'deb http://cran.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list && \
    apt-key add /jranke.key && \
    rm /jranke.key && \
    apt-get update && \
    apt-get install -f -y -t stretch-cran35 littler r-cran-littler r-base r-base-dev r-recommended && \
    apt-get install -f -y gdebi-core && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# register Java 8 with R
RUN R CMD javareconf

# make subdirectory mda_install
# COPY is always done as root!!!!
# so chmod for all
RUN mkdir /mda_install
COPY installations /mda_install
RUN chmod -R a+r /mda_install

# perform R installs for packages
# Separately, so if GitHub use runs out, restart has existing layers
RUN Rscript /mda_install/installNGCHM.R
RUN Rscript /mda_install/installForMBatch0.R
RUN Rscript /mda_install/installForMBatch1.R
RUN Rscript /mda_install/installForMBatch2.R
RUN Rscript /mda_install/installForMBatch3.R

RUN java -version && \
    R --version && \
    Rscript -e "installed.packages()"

## build with: docker build -t debian_r_java .
## rename/tag for push with: docker tag debian_r_java mdabcb/debian_r_java:2020-09-11-1000
## login to docker hub with: docker login --username=<your-username>
## repo may need be created on Docker Hub before pushing (depends on versions of stuff)
## push to docker hub with: docker push mdabcb/debian_r_java:2020-09-11-1000

