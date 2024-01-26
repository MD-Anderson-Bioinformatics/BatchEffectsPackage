FROM fedora:38

# docker build -t mdabcb/bcb_base:2024-01-11-0900 .
# Login to docker hub with: docker login
# repo may need be created on Docker Hub before pushing (depends on versions of stuff)
# Push to Docker Hub with: docker push mdabcb/bcb_base:2024-01-11-0900

# reminder, use {} around environmental variables, otherwise docker uses it as a literal

LABEL edu.mda.bcb.name="BCBBase" \
      edu.mda.bcb.sub="coj" \
      edu.mda.bcb.coj.version="2024-01-11-0900" \
      edu.mda.bcb.coj.Rversion="4.x" \
      edu.mda.bcb.coj.Javaversion="17" \
      edu.mda.bcb.coj.Linuxversion="Fedora 38"

####
#### generic setup for OS
####

# set timezone to prevent R time warnings
RUN ln -snf /usr/share/zoneinfo/US/Central /etc/localtime && echo "US/Central" > /etc/timezone
ENV TZ=US/Central

####
#### Make group to allow later users to access files,
#### without duplicating files/dirs in layer with chmod/chown
####

RUN groupadd -g 131313 bcb_base_group

####
#### Install/Setup Java 17 with variables
####

RUN dnf upgrade -y && \
    dnf -y install java-17-openjdk java-17-openjdk-devel java-17-openjdk-headless && \
    update-alternatives --install /usr/bin/java java /usr/lib/jvm/java-17-openjdk/bin/java 1 && \
    update-alternatives --install /usr/bin/javac javac /usr/lib/jvm/java-17-openjdk/bin/javac 1 && \
    update-alternatives --set java /usr/lib/jvm/java-17-openjdk/bin/java && \
    update-alternatives --set javac /usr/lib/jvm/java-17-openjdk/bin/javac && \
    dnf upgrade -y && \
    dnf clean all

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
ENV LD_LIBRARY_PATH=/usr/lib64/:/usr/lib/jvm/java-17-openjdk/lib/server

####
#### Install/Setup R 4.x+ plus OS packages for MBatch
####

# installs for R-related elements (like cairo)
# and build elements like Git and Maven
# set enables repo names need to be all lowercase, despite some online docs indicating otherwise
RUN dnf upgrade -y && \
    dnf -y install dnf-plugins-core && \
    dnf -y install libstdc++ ffmpeg-free ffmpeg-free-devel unzip diffutils openssl-devel libxml2-devel cairo-devel libXt-devel udunits2-devel proj-devel geos-devel gdal-devel sqlite sqlite-devel perl-Tk git git-gui ant maven nodejs file wget harfbuzz fribidi harfbuzz-devel fribidi-devel libcurl-devel freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel tesseract tesseract-devel leptonica-devel cargo poppler-cpp-devel ImageMagick-c++-devel libwebp-devel librsvg2-devel libgit2-devel gsl gsl-devel && \
    dnf upgrade -y && \
    dnf clean all
# npm install -g express-generator -- check if needed, has unpatched critical vulnerability

# Install all TexLive
RUN dnf upgrade -y && \
    dnf install -y texlive-2022 texlive-framed && \
    dnf upgrade -y && \
    dnf clean all

# install pandoc after texlive, since it depends on texlive
RUN dnf upgrade -y && \
    dnf -y install pandoc && \
    dnf upgrade -y && \
    dnf clean all

# install R 4.x+
RUN dnf upgrade -y && \
    dnf install -y R-java-devel && \
    dnf upgrade -y && \
    dnf clean all

RUN echo 'LC_CTYPE="C"' >> /usr/lib64/R/etc/Renviron && \
    echo 'LC_TIME="C"' >> /usr/lib64/R/etc/Renviron && \
    echo 'LC_MESSAGES="C"' >> /usr/lib64/R/etc/Renviron && \
    echo 'LC_MONETARY="C"' >> /usr/lib64/R/etc/Renviron && \
    echo 'LC_PAPER="C"' >> /usr/lib64/R/etc/Renviron && \
    echo 'LC_MEASUREMENT="C"' >> /usr/lib64/R/etc/Renviron

# register Java with R
RUN R CMD javareconf

####
#### Install/Setup MBatch (root)
####

# copy installs
# COPY is always done as root!!!!
COPY installations /bcbsetup/.
# required for NGCHM-R install (throws bad flag problem but works)
RUN ln -s /usr/bin/hostid /bin/hostname
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

####
#### Base mapping dir, owned by special group
####

# add base /BEA dir and set ownership
RUN mkdir /BEA && \
    chown -R :bcb_base_group /BEA && \
    chmod -R 775 /BEA

####
#### Python installs
####

# bcbuser is not created yet, so make dir, and chown/cmod

RUN mkdir /home/bcbuser && \
    cd /home/bcbuser && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    mkdir /home/bcbuser/conda && \
    bash /home/bcbuser/Miniconda3-latest-Linux-x86_64.sh -b -p /home/bcbuser/conda -f && \
    . /home/bcbuser/conda/bin/activate && \
    conda init && \
    conda update -y conda && \
    conda init bash && \
    chown -R :bcb_base_group /home/bcbuser && \
    chmod -R 775 /home/bcbuser

# link to gendev environment
RUN . /home/bcbuser/conda/bin/activate && \
    conda init && \
    conda create -y -p /home/bcbuser/conda/envs/gendev && \
    conda activate /home/bcbuser/conda/envs/gendev && \
    ln -s /home/bcbuser/conda/envs/gendev /BEA/gendev

# LD_LIBRARY_PATH needed for Python to find libraries due to conda env issues
ENV LD_LIBRARY_PATH=/usr/lib64/:/home/bcbuser/conda/envs/gendev/lib:/usr/lib/jvm/java-17-openjdk/lib/server

# need chmod at end to support other users
RUN . /home/bcbuser/conda/bin/activate && \
    conda init && \
    conda activate /home/bcbuser/conda/envs/gendev && \
    conda install -y -c conda-forge python==3.11.* && \
    conda install -y -c conda-forge pyinstaller && \
    conda install -y -c conda-forge flask && \
    conda install -y -c conda-forge waitress && \
    conda install -y -c conda-forge pylint && \
    conda install -y -c conda-forge setuptools && \
    conda install -y -c conda-forge pandas && \
    conda install -y -c conda-forge numpy && \
    conda install -y -c conda-forge matplotlib && \
    conda install -y -c conda-forge scanpy && \
    conda install -y -c conda-forge pillow && \
    conda install -y -c conda-forge jsonpickle && \
    conda install -y -c conda-forge requests && \
    conda install -y -c conda-forge xmltodict && \
    conda install -y -c conda-forge cryptography && \
    conda install -y -c conda-forge urllib3 && \
    conda install -y -c conda-forge scipy && \
    conda clean --all --yes && \
    chmod -R 777 /home/bcbuser/conda/envs/gendev
#    pip install pex && \  -- check if needed, has unpatched high vulnerability

# add universal access to conda
RUN cat /etc/bashrc /bcbsetup/conda.txt >> /bcbsetup/bashrc && \
    cp /bcbsetup/bashrc /etc/bashrc

RUN echo 'RETICULATE_PYTHON="/home/bcbuser/conda/envs/gendev/bin/python3"' >> /usr/lib64/R/etc/Renviron

####
#### End of build steps
####

RUN R CMD config --all

RUN ls -l /bcbsetup && \
    java -version && \
    R --version && \
    Rscript -e "installed.packages()"

ENV PATH="/home/bcbuser/conda/bin:$PATH"

