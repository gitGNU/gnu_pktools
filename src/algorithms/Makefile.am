AM_LDFLAGS = $(GDAL_LDFLAGS) @AM_LDFLAGS@
AM_CXXFLAGS = -I$(top_srcdir)/src $(GDAL_CFLAGS) @AM_CXXFLAGS@

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
noinst_LIBRARIES = libalgorithms.a

# where to install the headers on the system
libalgorithms_adir = $(includedir)/algorithms

# the list of header files that belong to the library (to be installed later)
libalgorithms_a_HEADERS = Egcs.h  Filter2d.h  Filter.h  Histogram.h ConfusionMatrix.h myfann_cpp.h

# the sources to add to the library and to add to the source distribution
libalgorithms_a_SOURCES = $(libalgorithms_a_HEADERS) Egcs.cc  Filter2d.cc  Filter.cc  Histogram.cc ConfusionMatrix.cc
###############################################################################
