AM_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(GDAL_CFLAGS) @AM_CXXFLAGS@
AM_LDFLAGS = $(GSL_LIBS) $(GDAL_LDFLAGS) $(top_builddir)/src/algorithms/libalgorithms.la $(top_builddir)/src/imageclasses/libimageClasses.la $(top_builddir)/src/fileclasses/libfileClasses.la $(top_builddir)/src/base/libbase.la @AM_LDFLAGS@
LDADD = $(GSL_LIBS) $(GDAL_LDFLAGS) $(top_builddir)/src/algorithms/libalgorithms.la $(top_builddir)/src/imageclasses/libimageClasses.la $(top_builddir)/src/fileclasses/libfileClasses.la $(top_builddir)/src/base/libbase.la
###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build and install (the names of the final binaries)
bin_PROGRAMS = pkinfo pkcrop pkreclass pkdiff pkgetmask pksetmask pkcreatect pkdumpimg pkdumpogr pksieve pkstatascii pkstatogr pkegcs pkextract pkfillnodata pkfilter pkkalman pkfilterdem pkenhance pkfilterascii pkdsm2shadow pkcomposite pkndvi pkpolygonize pkascii2img pksvm pkfssvm pkascii2ogr pkeditogr

# the program to build but not install (the names of the final binaries)
#noinst_PROGRAMS =  pkxcorimg pkgeom

if USE_FANN
bin_PROGRAMS += pkann pkfsann pkregann
pkann_SOURCES = $(top_srcdir)/src/algorithms/myfann_cpp.h pkann.cc
pkann_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(FANN_CFLAGS) -I$(top_srcdir)/src/algorithms $(AM_CXXFLAGS)
pkann_LDADD = $(FANN_LIBS) $(FANN_CFLAGS) $(AM_LDFLAGS)
pkfsann_SOURCES = $(top_srcdir)/src/algorithms/myfann_cpp.h $(top_srcdir)/src/algorithms/CostFactory.h pkfsann.h pkfsann.cc
pkfsann_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(FANN_CFLAGS) -I$(top_srcdir)/src/algorithms $(AM_CXXFLAGS)
pkfsann_LDADD = $(GSL_LIBS) $(FANN_LIBS) $(FANN_CFLAGS) $(AM_LDFLAGS) -lalgorithms
pkregann_SOURCES = $(top_srcdir)/src/algorithms/myfann_cpp.h pkregann.cc
pkregann_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(FANN_CFLAGS) -I$(top_srcdir)/src/algorithms $(AM_CXXFLAGS)
pkregann_LDADD = $(GSL_LIBS) $(FANN_LIBS) $(FANN_CFLAGS) $(AM_LDFLAGS)
endif

if USE_LAS
bin_PROGRAMS += pklas2img
pklas2img_SOURCES = pklas2img.cc
pklas2img_LDADD = $(top_srcdir)/src/lasclasses/liblasClasses.la -llas $(AM_LDFLAGS)
endif
if USE_NLOPT
bin_PROGRAMS += pkoptsvm
pkoptsvm_SOURCES = $(top_srcdir)/src/algorithms/OptFactory.h pkoptsvm.cc
pkoptsvm_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lnlopt
endif
# list of sources for the binaries
pkinfo_SOURCES = pkinfo.cc
pkcrop_SOURCES = pkcrop.cc
pkreclass_SOURCES = pkreclass.cc
pkdiff_SOURCES = pkdiff.cc
pkgetmask_SOURCES = pkgetmask.cc
pksetmask_SOURCES = pksetmask.cc
pkcreatect_SOURCES = pkcreatect.cc
pkdumpimg_SOURCES = pkdumpimg.cc
pkdumpogr_SOURCES = pkdumpogr.h pkdumpogr.cc
pksieve_SOURCES = pksieve.cc
pkstatascii_SOURCES = $(top_srcdir)/src/algorithms/StatFactory.h pkstatascii.cc
pkstatascii_LDADD = $(GSL_LIBS) $(AM_LDFLAGS)
pkstatogr_SOURCES = pkstatogr.cc
pkegcs_SOURCES = pkegcs.cc
pkegcs_LDADD = -lgdal $(AM_LDFLAGS) -lgdal
pkextract_SOURCES = pkextract.cc
pkfillnodata_SOURCES = pkfillnodata.cc
pkfilter_SOURCES = pkfilter.cc
pkfilter_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lgsl -lgdal
pkkalman_SOURCES = pkkalman.cc $(top_srcdir)/src/algorithms/ImgRegression.h $(top_srcdir)/src/algorithms/ImgRegression.cc
pkkalman_LDADD = -lalgorithms $(GSL_LIBS) $(AM_LDFLAGS) -lalgorithms
pkfilterdem_SOURCES = pkfilterdem.cc
pkenhance_SOURCES = pkenhance.cc
pkenhance_LDADD = $(AM_LDFLAGS) -lgdal
pkfilterascii_SOURCES = pkfilterascii.cc
pkfilterascii_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lgsl -lgdal
pkdsm2shadow_SOURCES = pkdsm2shadow.cc
pkcomposite_SOURCES = pkcomposite.cc
pkndvi_SOURCES = pkndvi.cc
pkpolygonize_SOURCES = pkpolygonize.cc
pksvm_SOURCES = $(top_srcdir)/src/algorithms/svm.h $(top_srcdir)/src/algorithms/svm.cpp pksvm.cc
pkfssvm_SOURCES = $(top_srcdir)/src/algorithms/svm.h $(top_srcdir)/src/algorithms/FeatureSelector.h  $(top_srcdir)/src/algorithms/CostFactory.h $(top_srcdir)/src/algorithms/svm.cpp pkfssvm.h pkfssvm.cc
pkfssvm_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lalgorithms
pkascii2img_SOURCES = pkascii2img.cc
pkascii2ogr_SOURCES = pkascii2ogr.cc
pkeditogr_SOURCES = pkeditogr.cc
pkeditogr_LDADD = $(AM_LDFLAGS)
###############################################################################
