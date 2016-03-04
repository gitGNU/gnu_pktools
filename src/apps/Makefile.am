AM_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(GDAL_CFLAGS) @AM_CXXFLAGS@
AM_LDFLAGS = $(GSL_LIBS) $(GDAL_LDFLAGS) $(top_builddir)/src/algorithms/libalgorithms.la $(top_builddir)/src/imageclasses/libimageClasses.la $(top_builddir)/src/fileclasses/libfileClasses.la @AM_LDFLAGS@
LDADD = $(GSL_LIBS) $(GDAL_LDFLAGS) $(top_builddir)/src/algorithms/libalgorithms.la $(top_builddir)/src/imageclasses/libimageClasses.la $(top_builddir)/src/fileclasses/libfileClasses.la
###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build and install (the names of the final binaries)
#bin_PROGRAMS = pkinfo pkcrop pkdiff pkgetmask pksetmask pkcreatect pkdumpimg pkdumpogr pksieve pkstat pkstatascii pkstatogr pkstatprofile pkegcs pkextract pkfillnodata pkfilter pkfilterdem pkfilterascii pkdsm2shadow pkcomposite pkpolygonize pkascii2img pksvm pkfssvm pkascii2ogr pkreclass
bin_PROGRAMS = pkinfo pkcrop pkdiff pkgetmask pksetmask pkcreatect pkdumpimg pkdumpogr pksieve pkstat pkstatascii pkstatogr pkstatprofile pkegcs pkextractimg pkextractogr pkfillnodata pkfilter pkfilterdem pkfilterascii pkdsm2shadow pkcomposite pkpolygonize pkascii2img pksvm pkfssvm pkascii2ogr pkreclass pkkalman

# the program to build but not install (the names of the final binaries)
#noinst_PROGRAMS = pkkalman


if USE_FANN
bin_PROGRAMS += pkann pkfsann pkregann
pkann_SOURCES = $(top_srcdir)/src/algorithms/myfann_cpp.h pkann.cc
pkann_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(FANN_CFLAGS) -I$(top_srcdir)/src/algorithms $(AM_CXXFLAGS)
pkann_LDADD = $(FANN_LIBS) $(FANN_CFLAGS) $(AM_LDFLAGS) -lgsl
pkfsann_SOURCES = $(top_srcdir)/src/algorithms/myfann_cpp.h $(top_srcdir)/src/algorithms/CostFactory.h pkfsann.h pkfsann.cc
pkfsann_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(FANN_CFLAGS) -I$(top_srcdir)/src/algorithms $(AM_CXXFLAGS)
pkfsann_LDADD = $(GSL_LIBS) $(FANN_LIBS) $(FANN_CFLAGS) $(AM_LDFLAGS) -lalgorithms -lgsl
pkregann_SOURCES = $(top_srcdir)/src/algorithms/myfann_cpp.h pkregann.cc
pkregann_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/base $(FANN_CFLAGS) -I$(top_srcdir)/src/algorithms $(AM_CXXFLAGS)
pkregann_LDADD = $(GSL_LIBS) $(FANN_LIBS) $(FANN_CFLAGS) $(AM_LDFLAGS)
endif

if USE_LAS
bin_PROGRAMS += pklas2img
pklas2img_SOURCES = pklas2img.cc
pklas2img_LDADD = -lgsl $(top_srcdir)/src/lasclasses/liblasClasses.la -llas $(AM_LDFLAGS)
endif
if USE_NLOPT
bin_PROGRAMS += pkoptsvm
pkoptsvm_SOURCES = $(top_srcdir)/src/algorithms/OptFactory.h $(top_srcdir)/src/algorithms/CostFactorySVM.h pkoptsvm.cc
pkoptsvm_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lnlopt -lgsl
endif
# list of sources for the binaries
pkinfo_SOURCES = pkinfo.cc
pkinfo_LDADD = $(AM_LDFLAGS) -lgsl
pkcrop_SOURCES = pkcrop.cc
pkcrop_LDADD = -lgsl
pkdiff_SOURCES = pkdiff.cc
pkdiff_LDADD = -lgsl
pkgetmask_SOURCES = pkgetmask.cc
pkgetmask_LDADD = $(GSL_LIBS)
pksetmask_SOURCES = pksetmask.cc
pksetmask_LDADD = $(GSL_LIBS)
pkcreatect_SOURCES = pkcreatect.cc
pkdumpimg_SOURCES = pkdumpimg.cc
pkdumpogr_SOURCES = pkdumpogr.h pkdumpogr.cc
pksieve_SOURCES = pksieve.cc
pkstat_SOURCES = $(top_srcdir)/src/algorithms/StatFactory.h $(top_srcdir)/src/algorithms/ImgRegression.h $(top_srcdir)/src/algorithms/ImgRegression.cc pkstat.cc
pkstat_LDADD = -lgsl -lgslcblas $(GSL_LIBS) $(AM_LDFLAGS)
pkstatprofile_SOURCES = $(top_srcdir)/src/algorithms/StatFactory.h pkstatprofile.cc
pkstatprofile_LDADD = -lgsl -lgslcblas $(GSL_LIBS) $(AM_LDFLAGS)
pkstatascii_SOURCES = $(top_srcdir)/src/algorithms/StatFactory.h pkstatascii.cc
pkstatascii_LDADD = -lgsl -lgslcblas $(GSL_LIBS) $(AM_LDFLAGS)
pkstatogr_SOURCES = pkstatogr.cc
pkstatogr_LDADD = -lgsl
pkegcs_SOURCES = pkegcs.cc
pkegcs_LDADD = $(AM_LDFLAGS) -lgsl -lgdal
pkextractimg_SOURCES = pkextractimg.cc
pkextractimg_LDADD = $(AM_LDFLAGS) -lgsl -lgdal
pkextractogr_SOURCES = pkextractogr.cc
pkextractogr_LDADD = $(AM_LDFLAGS) -lgsl -lgdal
pkfillnodata_SOURCES = pkfillnodata.cc
pkfilter_SOURCES = pkfilter.cc
pkfilter_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lgsl -lgdal
pkkalman_SOURCES = pkkalman.cc $(top_srcdir)/src/algorithms/ImgRegression.h $(top_srcdir)/src/algorithms/ImgRegression.cc
pkkalman_LDADD = -lalgorithms -lgsl $(GSL_LIBS) $(AM_LDFLAGS) -lalgorithms
pkfilterdem_SOURCES = pkfilterdem.cc
pkfilterdem_LDADD = -lgsl
pkfilterascii_SOURCES = pkfilterascii.cc
pkfilterascii_LDADD = $(GSL_LIBS) $(AM_LDFLAGS) -lgsl -lgdal
pkdsm2shadow_SOURCES = pkdsm2shadow.cc
pkdsm2shadow_LDADD = -lgsl
pkcomposite_SOURCES = pkcomposite.cc
pkpolygonize_SOURCES = pkpolygonize.cc
pksvm_SOURCES = $(top_srcdir)/src/algorithms/svm.h $(top_srcdir)/src/algorithms/svm.cpp pksvm.cc
pksvm_LDADD = -lgsl
pkfssvm_SOURCES = $(top_srcdir)/src/algorithms/svm.h $(top_srcdir)/src/algorithms/FeatureSelector.h  $(top_srcdir)/src/algorithms/CostFactorySVM.h $(top_srcdir)/src/algorithms/svm.cpp pkfssvm.cc
pkfssvm_LDADD = -lgsl $(GSL_LIBS) $(AM_LDFLAGS) -lalgorithms
pkascii2img_SOURCES = pkascii2img.cc
pkascii2ogr_SOURCES = pkascii2ogr.cc
pkreclass_SOURCES = pkreclass.cc
###############################################################################
