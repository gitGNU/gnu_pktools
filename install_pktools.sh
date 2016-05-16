########################################################################################
# install_pktools.sh script to install pktools on a debian based Linux distro (e.g., Ubuntu)
# Copyright (C) 2014 Pieter Kempeneers
#
# Basic usage: bash install_pktools.sh
# Advanced usage: bash install_pktools.sh [--enable-las] [--enable-nlopt] [--enable-fann]
#
# pktools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this script.  If not, see <http://www.gnu.org/licenses/>.
########################################################################################

#Set following variables to 0 if you do not want support for
LAS=1 #pklas2img
FANN=1 #pkann and pkfsann
NLOPT=1 #pkoptsvm

CONFIGURE="cmake .."

#!/bin/bash
function usage
{
    echo "usage: $0 [--enable-las] [--enable-fann] [-h]"
}

while [ "$1" != "" ]; do
    case $1 in
        --enable-las )             LAS=1
	    echo "las enabled"
                                ;;
        --enable-fann )             FANN=1
	    echo "fann enabled"
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

apt-get -y install software-properties-common
#get up to date
apt-get update


PPA_GDAL=$(apt-cache search 'libgdal-dev')

if [ -z "$PPA_GDAL" ];then
    add-apt-repository -y ppa:ubuntugis/ubuntugis-unstable
    #get up to date with new repository
    apt-get update
fi
    
#Install required pre-requisites for pktools
apt-get -y install g++ make gdal-bin libgdal-dev libgsl0-dev libarmadillo-dev

if [ "${FANN}" -eq 1 ];then
    #Install optional pre-requisites for Artificial Neural Network support
    apt-get -y install libfann-dev
    CONFIGURE="$CONFIGURE -DBUILD_WITH_FANN=ON"
fi

if [ "${LAS}" -eq 1 ];then
    #We do a manual installation of liblas to allow support for LAZ

    #install LASZIP
    cd /tmp
    wget https://github.com/LASzip/LASzip/releases/download/v2.2.0/laszip-src-2.2.0.tar.gz
    tar xzvf laszip-src-2.2.0.tar.gz
    cd laszip*
    ./configure
    make
    make install
    rm -rf /tmp/laszip*

    #hack: include files of optional package belong in their own directory
    mkdir -p /usr/local/include/laszip
    cd /usr/local/include/laszip
    for file in ../laszip*.hpp ../lasunzip*.hpp;do sudo ln -s $file $(basename $file);done

    #install liblas

    PPA_GTIFF=$(apt-cache search 'libgeotiff-dev')

    if [ -z "$PPA_GTIFF" ];then
	add-apt-repository -y ppa:ubuntugis/ubuntugis-unstable
	#get up to date with new repository
	apt-get update
    fi

    #The following packages need to be installed to support LASZIP within liblas 
   apt-get install -y cmake libgeotiff-dev libboost-program-options-dev libboost-dev libboost-thread-dev libboost-iostreams-dev libboost-filesystem-dev

    #get liblas 1.8.0 and install using cmake
    cd /tmp
    wget http://download.osgeo.org/liblas/libLAS-1.8.0.tar.bz2
    tar xjvf libLAS-1.8.0.tar.bz2
    cd libLAS-1.8.0
    mkdir makefiles
    cd makefiles
    cmake -G "Unix Makefiles" \
	-D GEOTIFF_INCLUDE_DIR=/usr/include/geotiff \
	-D LASZIP_INCLUDE_DIR=/usr/local/include \
	-D TIFF_INCLUDE_DIR=/usr/include \
	-D WITH_GEOTIFF=ON \
	-D WITH_LASZIP=ON \
	../
    make
    make install
    cd /tmp
    rm -rf /tmp/libLAS*
    #Install optional pre-requisites for libLAS support
    #apt-get -y install libboost-dev liblas-dev liblas-c-dev liblas1 liblas2 liblas-c2 python-liblas
    CONFIGURE="$CONFIGURE -DBUILD_WITH_LIBLAS=ON"
fi

ldconfig

#Install pktools
cd /tmp
#todo: use package manager once pktools in repository
#Manual installation
wget --progress=dot:mega "download.savannah.gnu.org/releases/pktools/pktools-latest.tar.gz"
tar xzvf pktools-latest.tar.gz
cd pktools-*
mkdir build
cd build
$CONFIGURE
#./configure --enable-fann --enable-las
make
make install
ldconfig

rm -rf /tmp/pktools-*
