# R and Java for MDA Batch Effects Projects

This image includes AlmaLinux 8, R 4+, and Java 8 along with other required packages for use with MD Anderson Cancer Center Bioinformatics and Computational Biology's Batch Effects effort.

This is for educational and research purposes only. 

Additional information on Batch Effects can be found at http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview
Downloads and details on Standardized Data are available at http://bioinformatics.mdanderson.org/TCGA/databrowser/

This can be built with: 
docker build -t bcb_base .

Rename/tag for push with: 
docker tag bcb_base mdabcb/bcb_base:2021-08-06-1100

Login to docker hub with: 
docker login

Push to Docker Hub with: 
docker push mdabcb/bcb_base:2021-08-06-1100

