AM_LDFLAGS = $(GDAL_LDFLAGS) @AM_LDFLAGS@ -llas
AM_CXXFLAGS = -I$(top_srcdir)/src $(GDAL_CFLAGS) @AM_CXXFLAGS@

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
#noinst_LIBRARIES = libfileClasses.a
lib_LTLIBRARIES = libfileClasses.la

# where to install the headers on the system
libfileClasses_ladir = $(includedir)/fileclasses

# the list of header files that belong to the library (to be installed later)
libfileClasses_la_HEADERS = FileReaderAscii.h

# the sources to add to the library and to add to the source distribution
libfileClasses_la_SOURCES = $(libfileClasses_la_HEADERS) FileReaderAscii.cc

if USE_LAS
# the list of header files that belong to the library (to be installed later)
libfileClasses_la_HEADERS += FileReaderLas.h
# the sources to add to the library and to add to the source distribution
libfileClasses_la_SOURCES += FileReaderLas.h FileReaderLas.cc
endif

###############################################################################
