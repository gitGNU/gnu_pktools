AM_LDFLAGS = $(GSL_LIBS) $(GDAL_LDFLAGS) @AM_LDFLAGS@ 
AM_CXXFLAGS = -I$(top_srcdir)/src $(GDAL_CFLAGS) @AM_CXXFLAGS@

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
#do not want to install pktestoption
noinst_PROGRAMS = pktestProspect

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
noinst_LIBRARIES = libmodels.a

# where to install the headers on the system
libmodels_adir = $(includedir)/models

# the list of header files that belong to the library (to be installed later)
libmodels_a_HEADERS = Prospect.h

if USE_NLOPT
libmodels_a_HEADERS += SensorModel.h
endif

# the sources to add to the library and to add to the source distribution
libmodels_a_SOURCES = $(libmodels_a_HEADERS) Prospect.cc 
###############################################################################

# list of sources for the binaries
pktestProspect_SOURCES = pktestProspect.cc
pktestProspect_LDADD = $(top_builddir)/src/models/libmodels.a $(top_builddir)/src/models/libprospect.a libf95.a
