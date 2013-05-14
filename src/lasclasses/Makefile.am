AM_LDFLAGS = $(GDAL_LDFLAGS) @AM_LDFLAGS@ -llas
AM_CXXFLAGS = -I$(top_srcdir)/src $(GDAL_CFLAGS) @AM_CXXFLAGS@

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
lib_LTLIBRARIES = liblasClasses.la

# where to install the headers on the system
liblasClasses_ladir = $(includedir)/lasclasses

# the list of header files that belong to the library (to be installed later)
liblasClasses_la_HEADERS = FileReaderLas.h

# the sources to add to the library and to add to the source distribution
liblasClasses_la_SOURCES = $(liblasClasses_la_HEADERS) FileReaderLas.cc

###############################################################################
