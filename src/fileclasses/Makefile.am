AM_LDFLAGS = $(GDAL_LDFLAGS) @AM_LDFLAGS@
AM_CXXFLAGS = -I$(top_srcdir)/src $(GDAL_CFLAGS) @AM_CXXFLAGS@

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
lib_LTLIBRARIES = libfileClasses.la

# where to install the headers on the system
libfileClasses_ladir = $(includedir)/pktools/fileclasses

## Instruct libtool to include ABI version information in the generated shared
## library file (.so).  The library ABI version is defined in configure.ac, so
## that all version information is kept in one place.
libfileClasses_la_LDFLAGS = -version-info $(PKTOOLS_SO_VERSION) $(AM_LDFLAGS)

# the list of header files that belong to the library (to be installed later)
libfileClasses_la_HEADERS = FileReaderAscii.h

# the sources to add to the library and to add to the source distribution
libfileClasses_la_SOURCES = $(libfileClasses_la_HEADERS) FileReaderAscii.cc

###############################################################################
