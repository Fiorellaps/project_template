#!/bin/bash

#Creacion del directorio para las librerias

mkdir software software/deps_r
chmod 755 software/deps_r
chmod 777 -R  software/deps_r/
Rscript -e "Sys.setenv(R_INSTALL_STAGED = FALSE)"

#Instalacion de paquetes necesarios para linux

sudo apt-get install curl libssl-dev libcurl4-openssl-dev libxml2-dev

#Instalacion de librerias R

PKGS_CRAN="'readxl' 'igraph' 'linkcomm' 'data.table' 'BiocManager' 'devtools' 'ggplot2' 'zoo' 'dplyr' 'data.table' 'pkgload' 'xml2'"
PKGS_BM="'STRINGdb' 'clusterProfiler' 'DOSE' 'biomaRt' 'enrichplot' 'org.Hs.eg.db'"

for PKG_CRAN in $PKGS_CRAN
    do
        Rscript -e 'install.packages('$PKG_CRAN', repos="https://cran.rstudio.com/", lib="software/deps_r" )'
        Rscript -e 'library('$PKG_CRAN', lib.loc="software/deps_r")'
    done

Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.12",ask=FALSE)'

for PKG_BM in $PKGS_BM
    do
        Rscript -e 'BiocManager::install('$PKG_BM', lib="software/deps_r" )'
        Rscript -e 'library('$PKG_BM',lib.loc="software/deps_r")'
    done
    
#A�adir el path a la variable libPaths

Rscript -e '.libPaths(c(.libPaths(),paste(getwd(),"software/deps_r",sep="/")))'
echo "************************** ALL DONE ********************************"
echo "Librerias cargadas en el directorio /software/deps_r"