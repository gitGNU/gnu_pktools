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

# the library names to build (note we are building static libs only)
noinst_LIBRARIES = libalgorithms.a

# where to install the headers on the system
libalgorithms_adir = $(includedir)/algorithms

# the list of header files that belong to the library (to be installed later)
libalgorithms_a_HEADERS = Egcs.h Filter2d.h Filter.h StatFactory.h ConfusionMatrix.h svm.h FeatureSelector.h

if USE_FANN
libalgorithms_a_HEADERS += myfann_cpp.h
endif

if USE_NLOPT
libalgorithms_a_HEADERS += SensorModel.h OptFactory.h
endif

# the sources to add to the library and to add to the source distribution
libalgorithms_a_SOURCES = $(libalgorithms_a_HEADERS) Egcs.cc Filter2d.cc Filter.cc ConfusionMatrix.cc svm.cpp
###############################################################################

# list of sources for the binaries
#pktestStat_SOURCES = pktestStat.cc
#pktestStat_LDADD = $(GSL_LIBS) -lgslwrap $(top_builddir)/src/algorithms/libalgorithms.a $(top_builddir)/src/imageclasses/libimageClasses.a -lgslwrap
