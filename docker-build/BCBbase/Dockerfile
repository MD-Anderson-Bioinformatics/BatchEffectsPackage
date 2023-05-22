FROM almalinux:8

# cd /code/BatchEffects_clean/BatchEffectsPackage/docker-build/BCBbase/
# docker build -t mdabcb/bcb_base:2023-03-09-1313 .
# Login to docker hub with: docker login
# repo may need be created on Docker Hub before pushing (depends on versions of stuff)
# Push to Docker Hub with: docker push mdabcb/bcb_base:2023-03-09-1313

# save
# docker save -o ~/bcb_base_22023-03-09-1313.tar mdabcb/bcb_base:2023-03-09-1313

# reminder, use {} around environmental variables, otherwise docker uses it as a literal

LABEL edu.mda.bcb.name="BCBBase" \
      edu.mda.bcb.sub="coj" \
      edu.mda.bcb.coj.version="2023-03-09-1313" \
      edu.mda.bcb.coj.Rversion="4.x" \
      edu.mda.bcb.coj.Javaversion="17" \
      edu.mda.bcb.coj.AlmaLinuxversion="8"


####
#### generic setup for OS
####

# set timezone to prevent R time warnings
RUN ln -snf /usr/share/zoneinfo/US/Central /etc/localtime && echo "US/Central" > /etc/timezone
ENV TZ=US/Central

####
#### Install/Setup Java 17 with variables
####

RUN dnf -y install java-17-openjdk java-17-openjdk-devel java-17-openjdk-headless && \
    update-alternatives --install /usr/bin/java java /usr/lib/jvm/java-17-openjdk/bin/java 1 && \
    update-alternatives --install /usr/bin/javac javac /usr/lib/jvm/java-17-openjdk/bin/javac 1 && \
    update-alternatives --set java /usr/lib/jvm/java-17-openjdk/bin/java && \
    update-alternatives --set javac /usr/lib/jvm/java-17-openjdk/bin/javac && \
    yum update -y && \
    yum clean all

# adding these variables fixes issue where registering Java with R (keeps getting "lost")
ENV JAR=/usr/lib/jvm/java-17-openjdk/bin/jar 
ENV JAVAH=/usr/lib/jvm/java-17-openjdk/bin/javac 
ENV JAVA_LIBS="-L/usr/lib/jvm/java-17-openjdk/lib/server -ljvm" 
ENV JAVA_CPPFLAGS="-I/usr/lib/jvm/java-17-openjdk/include -I/usr/lib/jvm/java-17-openjdk/include/linux" 
ENV JAVA_LD_LIBRARY_PATH=/usr/lib/jvm/java-17-openjdk/lib/server
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk/
ENV PATH=/usr/lib/jvm/java-17-openjdk/bin:$PATH
ENV JAVA=/usr/lib/jvm/java-17-openjdk/bin/java
ENV JAVAC=/usr/lib/jvm/java-17-openjdk/bin/javac
ENV JAVAH=/usr/lib/jvm/java-17-openjdk/bin/javac
ENV JAR=/usr/lib/jvm/java-17-openjdk/bin/jar
ENV LD_LIBRARY_PATH=/usr/lib/jvm/java-17-openjdk/lib/server

####
#### Install/Setup R 4.x+ plus OS packages for MBatch
####

# installs for R-related elements (like cairo)
# and build elements like Git and Maven
# set enables repo names need to be all lowercase, despite some online docs indicating otherwise
RUN dnf -y install dnf-plugins-core && \
    dnf config-manager --set-enabled powertools && \
    dnf config-manager --set-enabled appstream && \
    dnf -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    yum install -y https://download1.rpmfusion.org/free/el/rpmfusion-free-release-8.noarch.rpm https://download1.rpmfusion.org/nonfree/el/rpmfusion-nonfree-release-8.noarch.rpm && \
    yum install -y http://rpmfind.net/linux/centos/8-stream/PowerTools/x86_64/os/Packages/SDL2-2.0.10-2.el8.x86_64.rpm && \
    yum install -y ffmpeg ffmpeg-devel && \
    yum install -y unzip diffutils openssl-devel libxml2-devel cairo-devel libXt-devel && \
    yum install -y udunits2-devel proj-devel geos-devel gdal-devel sqlite sqlite-devel && \
    yum install -y perl-Tk git git-gui ant maven nodejs file wget && \
    yum install -y harfbuzz fribidi harfbuzz-devel fribidi-devel ffmpeg-devel && \
    yum install -y freetype-devel libpng-devel libtiff-devel && \
    yum install -y libcurl-devel libjpeg-turbo-devel && \
    yum install -y tesseract tesseract-devel leptonica-devel cargo && \
    yum install -y poppler-cpp-devel ImageMagick-c++-devel && \
    yum install -y libwebp-devel librsvg2-devel libgit2-devel && \
    yum install -y gsl gsl-devel && \
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

####
#### Requirements for MBatch Python install/setup
####

RUN mkdir -p /home/bcbuser && \
    cd /home/bcbuser && \
    wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh && \
    mkdir -p /home/bcbuser/conda && \
    bash /home/bcbuser/Anaconda3-2022.10-Linux-x86_64.sh -b -p /home/bcbuser/conda -f && \
    source /home/bcbuser/conda/bin/activate && \
    conda init && \
    conda update -y conda

# LD_LIBRARY_PATH needed for Python to find libraries due to conda env issues
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bcbuser/conda/envs/gendev/lib

RUN source /home/bcbuser/conda/bin/activate && \
    conda init && \
    conda create -y -p /home/bcbuser/conda/envs/gendev && \
    conda activate /home/bcbuser/conda/envs/gendev && \
    conda install -y -c conda-forge python==3.9 && \
    conda install -y -c conda-forge tensorflow && \
    conda install -y -c conda-forge pandas && \
    conda install -y -c conda-forge numpy && \
    conda install -y -c conda-forge scipy && \
    conda install -y -c conda-forge pyinstaller && \
    conda install -y -c conda-forge mypy && \
    conda install -y -c conda-forge pylint && \
    conda install -y -c conda-forge flask && \
    conda install -y -c conda-forge waitress && \
    conda install -y -c conda-forge setuptools && \
    conda install -y -c conda-forge tox && \
    conda install -y -c conda-forge pipreqs && \
    conda install -y -c conda-forge matplotlib && \
    conda install -y -c conda-forge nptyping && \
    conda install -y -c conda-forge pillow && \
    conda install -y -c conda-forge jsonpickle && \
    conda install -y -c conda-forge xmltodict && \
    pip install pex

####
#### Requirements for MBatch R install/setup
####

# copy installs
# COPY is always done as root!!!!
COPY installations /bcbsetup/.
# perform R installs for packages
# Separately, so if GitHub use runs out, restart has existing layers
RUN Rscript /bcbsetup/installNGCHM.R
RUN Rscript /bcbsetup/installForMBatch0.R
RUN Rscript /bcbsetup/installForMBatch1.R
RUN Rscript /bcbsetup/installForMBatch2.R
RUN Rscript /bcbsetup/installForMBatch3.R
RUN Rscript /bcbsetup/postInstall.R
# needed to support other linkouts
RUN cp /bcbsetup/custom.js /usr/lib64/R/library/NGCHM/extdata/custom.js && \
    chmod 755 /usr/lib64/R/library/NGCHM/extdata/custom.js

# add universal access to conda
RUN cat /etc/bashrc /bcbsetup/conda.txt >> /bcbsetup/bashrc && \
    cp /bcbsetup/bashrc /etc/bashrc

####
#### End of build steps
####

RUN R CMD config --all

RUN ls -l /bcbsetup && \
    java -version && \
    R --version && \
    Rscript -e "installed.packages()"

ENV PATH="/home/bcbuser/conda/bin:$PATH"

