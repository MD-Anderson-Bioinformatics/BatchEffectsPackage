# Debian, R, and Java for MDA Batch Effects Projects

This image includes Debian 9, R 3.6+, and Java 8 along with other required packages for use with MD Anderson Cancer Center Bioinformatics and Computational Biology's Batch Effects effort.

This is for educational and research purposes only. 

Additional information on Batch Effects can be found at http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview
Downloads and details on Standardized Data are available at http://bioinformatics.mdanderson.org/TCGA/databrowser/


This can be built with: docker build -t debian_r_java .
Rename/tag for push with: docker tag debian_r_java mdabcb/debian_r_java:2019-06-24-1600
Login to docker hub with: docker login
Push to Docker Hub with: docker push mdabcb/debian_r_java:2019-06-24-1600

