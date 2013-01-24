# List of pktool applications (command line utilities)

<!-- To create the current list of pktool apps, you can use:
find src/apps/ -name 'pk*.cc'|sed -e 's/.*\//- /' -e 's/\.cc//'|sort 
-->

- pkascii2img
- pkascii2ogr
- pkclassify_nn
- pkclassify_svm
- pkcreatect
- pkcrop
- pkdiff
- pkdsm2shadow
- pkdumpimg
- pkdumpogr
- pkegcs
- pkextract
- pkfillnodata
- pkfilter
- pkfs_nn
- pkfs_svm
- pkgeom
- pkgetchandelier
- pkgetmask
- pkinfo
- pklas2img
- pkmosaic
- pkndvi
- pkopt_svm
- pkpolygonize
- pkreclass
- pksensormodel
- pksetchandelier
- pksetmask
- pksieve
- pkstat
- pkstatogr
- pkxcorimg

/*! \page pkinfo pkinfo

lists specific information about a raster dataset (similar to gdalinfo but information given is driven by options)

##pkinfo_synopsis SYNOPSIS

pkinfo [--help-general] [-mm] [-stats] [-hist] [-nogcp] [-nomd]
         [-noct] [-nofl] [-checksum] [-proj4] [-mdd domain]*
	 [-sd subdataset] datasetname

##pkinfo_description DESCRIPTION
