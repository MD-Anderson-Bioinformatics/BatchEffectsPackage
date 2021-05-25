# Debian, R, and Java for MDA Batch Effects Projects

This image includes CentOS 8, R 4+, and Java 8 along with other required packages for use with MD Anderson Cancer Center Bioinformatics and Computational Biology's Batch Effects effort.

This is for educational and research purposes only. 

Additional information on Batch Effects can be found at http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview
Downloads and details on Standardized Data are available at http://bioinformatics.mdanderson.org/TCGA/databrowser/

This can be built with: docker build -t centos_r_java .
Rename/tag for push with: docker tag centos_r_java mdabcb/centos_r_java:2021-04-27-1030
Login to docker hub with: docker login
Push to Docker Hub with: docker push mdabcb/centos_r_java:2021-04-27-1030

