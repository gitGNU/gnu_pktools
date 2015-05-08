#!/bin/bash
pkcrop -i data/lena.tif -o data/output/lena_10.tif --ulx 256 -uly 256 -lrx 512 -lry 0
pkcrop -i data/lena.tif -o data/output/lena_00.tif --ulx 0 -uly 256 -lrx 256 -lry 0
pkcrop -i data/lena.tif -o data/output/lena_11.tif --ulx 256 -uly 512 -lrx 512 -lry 256
pkcrop -i data/lena.tif -o data/output/lena_01.tif --ulx 0 -uly 512 -lrx 256 -lry 256
pkcomposite -i data/output/lena_00.tif -i data/output/lena_01.tif -i data/output/lena_10.tif -i data/output/lena_11.tif -o data/output/lena.tif
pkdiff -ref data/lena.tif -i data/output/lena.tif

