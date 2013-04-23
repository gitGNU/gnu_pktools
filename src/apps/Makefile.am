AM_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(GDAL_CFLAGS) @AM_CXXFLAGS@
AM_LDFLAGS = $(GDAL_LDFLAGS) $(top_builddir)/src/algorithms/libalgorithms.la $(top_builddir)/src/imageclasses/libimageClasses.la $(top_builddir)/src/fileclasses/libfileClasses.la @AM_LDFLAGS@
#LDADD = $(GDAL_LDFLAGS) $(top_builddir)/src/algorithms/libalgorithms.la $(top_builddir)/src/models/libmodels.la $(top_builddir)/src/models/libprospect.la $(top_builddir)/src/imageclasses/libimageClasses.la
LDADD = $(GDAL_LDFLAGS) $(top_builddir)/src/algorithms/libalgorithms.la $(top_builddir)/src/models/libmodels.la $(top_builddir)/src/imageclasses/libimageClasses.la $(top_builddir)/src/fileclasses/libfileClasses.la
###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build and install (the names of the final binaries)
bin_PROGRAMS = pkinfo pkcrop pkreclass pkgetmask pksetmask pkcreatect pkdumpimg pkdumpogr pksieve pkstat pkstatogr pkegcs pkextract pkfillnodata pkfilter pkdsm2shadow pkmosaic pkndvi pkpolygonize pkascii2img pkdiff pkclassify_svm pkfs_svm pkascii2ogr pkeditogr

# the program to build but not install (the names of the final binaries)
#noinst_PROGRAMS =  pkxcorimg pkgeom

if USE_FANN
bin_PROGRAMS += pkclassify_nn pkfs_nn
pkclassify_nn_SOURCES = $(top_srcdir)/src/algorithms/myfann_cpp.h pkclassify_nn.h pkclassify_nn.cc
pkclassify_nn_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(FANN_CFLAGS) -I$(top_srcdir)/src/algorithms $(AM_CXXFLAGS)
pkclassify_nn_LDADD = $(FANN_LIBS) $(FANN_CFLAGS) $(AM_LDFLAGS)
pkfs_nn_SOURCES = $(top_srcdir)/src/algorithms/myfann_cpp.h pkfs_nn.cc
pkfs_nn_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(FANN_CFLAGS) -I$(top_srcdir)/src/algorithms $(AM_CXXFLAGS)
pkfs_nn_LDADD = $(GSL_LIBS) $(FANN_LIBS) $(FANN_CFLAGS) $(AM_LDFLAGS)
endif

if USE_LAS
bin_PROGRAMS += pklas2img
pklas2img_SOURCES = pklas2img.cc
pklas2img_LDADD = -llas $(AM_LDFLAGS)
endif
if USE_NLOPT
bin_PROGRAMS += pksensormodel pkopt_svm
pksensormodel_SOURCES = $(top_srcdir)/src/models/SensorModel.h pksensormodel.h pksensormodel.cc
pksensormodel_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lnlopt -lm -lgslwrap
pkopt_svm_SOURCES = $(top_srcdir)/src/algorithms/OptFactory.h pkclassify_nn.h pkopt_svm.cc
pkopt_svm_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lnlopt
endif
# list of sources for the binaries
pkinfo_SOURCES = pkinfo.cc
pkcrop_SOURCES = pkcrop.cc
pkreclass_SOURCES = pkreclass.cc
pkgetmask_SOURCES = pkgetmask.cc
pksetmask_SOURCES = pksetmask.cc
pkcreatect_SOURCES = pkcreatect.cc
pkdumpimg_SOURCES = pkdumpimg.cc
pkdumpogr_SOURCES = pkdumpogr.h pkdumpogr.cc
pksieve_SOURCES = pksieve.cc
pkstat_SOURCES = $(top_srcdir)/src/algorithms/StatFactory.h pkstat.cc
pkstat_LDADD = -llas $(GSL_LIBS) $(AM_LDFLAGS)
pkstatogr_SOURCES = pkstatogr.cc
pkegcs_SOURCES = pkegcs.cc
pkegcs_LDADD = -lgdal $(AM_LDFLAGS) -lgdal
pkextract_SOURCES = pkextract.cc
pkfillnodata_SOURCES = pkfillnodata.cc
pkfilter_SOURCES = pkfilter.cc
pkfilter_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lgsl -lgslwrap -lgdal
pkdsm2shadow_SOURCES = pkdsm2shadow.cc
pkmosaic_SOURCES = pkmosaic.cc
pkndvi_SOURCES = pkndvi.cc
pkpolygonize_SOURCES = pkpolygonize.cc
pkdiff_SOURCES = pkdiff.cc
pkclassify_svm_SOURCES = $(top_srcdir)/src/algorithms/svm.h $(top_srcdir)/src/algorithms/svm.cpp pkclassify_nn.h pkclassify_svm.cc
pkfs_svm_SOURCES = $(top_srcdir)/src/algorithms/svm.h $(top_srcdir)/src/algorithms/svm.cpp pkclassify_nn.h pkfs_svm.cc
pkfs_svm_LDADD = $(GSL_LIBS) $(AM_LDFLAGS)
pkascii2img_SOURCES = pkascii2img.cc
pkascii2ogr_SOURCES = pkascii2ogr.cc
pkeditogr_SOURCES = pkeditogr.cc
#pkxcorimg_SOURCES= pkxcorimg.cc
#pkgeom_SOURCES = pkgeom.cc
###############################################################################
