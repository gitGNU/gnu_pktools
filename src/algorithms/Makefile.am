AM_LDFLAGS = $(GSL_LIBS) $(GDAL_LDFLAGS) @AM_LDFLAGS@ 
AM_CXXFLAGS = -I$(top_srcdir)/src $(GDAL_CFLAGS) @AM_CXXFLAGS@

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
#do not want to install pktestoption
#noinst_PROGRAMS = pktestStat

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

#noinst_LIBRARIES = libalgorithms.a
lib_LTLIBRARIES = libalgorithms.la

# where to install the headers on the system
libalgorithms_ladir = $(includedir)/pktools/algorithms

## Instruct libtool to include ABI version information in the generated shared
## library file (.so).  The library ABI version is defined in configure.ac, so
## that all version information is kept in one place.
libalgorithms_la_LDFLAGS = -version-info $(PKTOOLS_SO_VERSION) $(AM_LDFLAGS)

# the list of header files that belong to the library (to be installed later)
libalgorithms_la_HEADERS = Egcs.h Filter2d.h Filter.h StatFactory.h ConfusionMatrix.h svm.h CostFactory.h CostFactorySVM.h FeatureSelector.h

if USE_FANN
libalgorithms_la_HEADERS += myfann_cpp.h
endif

# if USE_NLOPT
# libalgorithms_la_HEADERS += OptFactory.h
# endif

# the sources to add to the library and to add to the source distribution
libalgorithms_la_SOURCES = $(libalgorithms_la_HEADERS) Egcs.cc Filter2d.cc Filter.cc ConfusionMatrix.cc svm.cpp CostFactorySVM.cc
###############################################################################

# list of sources for the binaries
#pktestStat_SOURCES = pktestStat.cc
#pktestStat_LDADD = $(GSL_LIBS) $(top_builddir)/src/algorithms/libalgorithms.la $(top_builddir)/src/imageclasses/libimageClasses.la
