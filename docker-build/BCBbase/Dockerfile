FROM almalinux:8

# BEA_VERSION_TIMESTAMP
## This can be built with: docker build -t mdabcb/bcb_base:BEA_VERSION_TIMESTAMP .
## Rename/tag for push with: docker tag bcb_base mdabcb/bcb_base:BEA_VERSION_TIMESTAMP
## Login to docker hub with: docker login
## repo may need be created on Docker Hub before pushing (depends on versions of stuff)
## Push to Docker Hub with: docker push mdabcb/bcb_base:BEA_VERSION_TIMESTAMP

# save
# docker save -o ~/bcb_base_BEA_VERSION_TIMESTAMP.tar mdabcb/bcb_base:BEA_VERSION_TIMESTAMP

# reminder, use {} around environmental variables, otherwise docker uses it as a literal

LABEL edu.mda.bcb.name="BCBBase" \
      edu.mda.bcb.sub="coj" \
      edu.mda.bcb.coj.version="BEA_VERSION_TIMESTAMP" \
      edu.mda.bcb.coj.Rversion="4.x" \
      edu.mda.bcb.coj.Javaversion="8" \
      edu.mda.bcb.coj.AlmaLinuxversion="8"

# Use Java 8, still supported, and all the tools depend on it -- simplifies other installs and maintenance
# also set environment variables needed by R and rJava, that are not otherwise preserved
# adding LD_LIBRARY_PATH fixes issue where registering Java with R (keeps getting "lost")
RUN dnf -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    dnf -y install java-1.8.0-openjdk java-1.8.0-openjdk-devel java-1.8.0-openjdk-headless && \
    update-alternatives --install /usr/bin/java java /usr/lib/jvm/java-1.8.0-openjdk/jre/bin/java 1 && \
    update-alternatives --install /usr/bin/javac javac /usr/lib/jvm/java-1.8.0-openjdk/bin/javac 1 && \
    update-alternatives --install /usr/bin/jar jar /usr/lib/jvm/java-1.8.0-openjdk/bin/jar 1 && \
    update-alternatives --install /usr/bin/javah javah /usr/lib/jvm/java-1.8.0-openjdk/bin/javah 1 && \
    echo "export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk/" >> /etc/profile.d/java.sh && \
    echo 'export PATH="/usr/lib/jvm/java-1.8.0-openjdk/bin:$PATH"' >> /etc/profile.d/java.sh && \
    echo "export JAVA=/usr/lib/jvm/java-1.8.0-openjdk/jre/bin/java" >> /etc/profile.d/rJava.sh && \
    echo "export JAVAC=/usr/lib/jvm/java-1.8.0-openjdk/bin/javac" >> /etc/profile.d/rJava.sh && \
    echo "export JAVAH=/usr/lib/jvm/java-1.8.0-openjdk/bin/javah" >> /etc/profile.d/rJava.sh && \
    echo "export JAR=/usr/lib/jvm/java-1.8.0-openjdk/bin/jar" >> /etc/profile.d/rJava.sh && \
    echo "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/lib/jvm/java-1.8.0-openjdk/lib/amd64:/usr/lib/jvm/java-1.8.0-openjdk/jre/lib/amd64/server" >> /etc/profile.d/r.sh && \
    source /etc/profile.d/java.sh && \
    source /etc/profile.d/rJava.sh && \
    source /etc/profile.d/r.sh && \
    yum update -y && \
    yum clean all

# installs for R-related elements (like cairo)
# and build elements like Git and Maven
# set enables repo names need to be all lowercase, despite some online docs indicating otherwise
RUN dnf -y install dnf-plugins-core && \
    dnf config-manager --set-enabled powertools && \
    dnf config-manager --set-enabled appstream && \
    yum install -y unzip diffutils openssl-devel libxml2-devel cairo-devel libXt-devel && \
    yum install -y udunits2-devel proj-devel geos-devel gdal-devel sqlite sqlite-devel && \
    yum install -y perl-Tk git git-gui ant maven nodejs file wget && \
    yum update -y && \
    yum clean all && \
    npm install -g express-generator

# Install all TexLive
RUN dnf install -y texliv*

# install pandoc after texlive, since it depends on texlive
RUN yum install -y pandoc && \
    yum update -y && \
    yum clean all

# install R 4.x+
RUN dnf install -y R-java-devel

RUN echo 'LC_CTYPE="C"' >> /usr/lib64/R/.Renviron && \
    echo 'LC_TIME="C"' >> /usr/lib64/R/.Renviron && \
    echo 'LC_MESSAGES="C"' >> /usr/lib64/R/.Renviron && \
    echo 'LC_MONETARY="C"' >> /usr/lib64/R/.Renviron && \
    echo 'LC_PAPER="C"' >> /usr/lib64/R/.Renviron && \
    echo 'LC_MEASUREMENT="C"' >> /usr/lib64/R/.Renviron

# register Java with R
RUN R CMD javareconf

# other required OS packages for R
RUN yum install -y libcurl-devel libjpeg-turbo-devel && \
    yum update -y && \
    yum clean all

# re-register Java with R (keeps getting "lost")
RUN source /etc/profile.d/java.sh && \
    source /etc/profile.d/rJava.sh && \
    source /etc/profile.d/r.sh && \
    R CMD javareconf

RUN yum install -y poppler-cpp-devel ImageMagick-c++-devel libwebp-devel librsvg2-devel libgit2-devel harfbuzz-devel fribidi-devel && \
    yum update -y && \
    yum clean all

# copy installs
# COPY is always done as root!!!!
COPY installations /bcbsetup/.
# perform R installs for packages
# Separately, so if GitHub use runs out, restart has existing layers
RUN Rscript /bcbsetup/installNGCHM.R
RUN Rscript /bcbsetup/installRJava.R
RUN Rscript /bcbsetup/installForMBatch0.R
RUN Rscript /bcbsetup/installForMBatch1.R
RUN Rscript /bcbsetup/installForMBatch2.R
RUN Rscript /bcbsetup/installForMBatch3.R
RUN Rscript /bcbsetup/postInstall.R

# needed to support other linkouts
RUN cp /bcbsetup/custom.js /usr/lib64/R/library/NGCHM/extdata/custom.js

RUN R CMD config --all

RUN ls -l /bcbsetup && \
    java -version && \
    R --version && \
    Rscript -e "installed.packages()"

