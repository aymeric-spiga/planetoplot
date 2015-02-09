#! /bin/bash

wget http://repo.continuum.io/miniconda/Miniconda-3.6.0-Linux-x86.sh 
chmod 755 Miniconda-3.6.0-Linux-x86.sh 
./Miniconda-3.6.0-Linux-x86.sh 
source $HOME/.bashrc
conda install numpy
conda install matplotlib
conda install netCDF4

##the following lib are not mandatory
#conda install basemap
#conda install pip
#conda install scipy
