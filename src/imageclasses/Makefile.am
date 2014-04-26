AM_LDFLAGS = $(GSL_LDFLAGS) $(GDAL_LDFLAGS) @AM_LDFLAGS@
AM_CXXFLAGS = -I$(top_srcdir)/src $(GSL_CFLAGS) $(GDAL_CFLAGS) @AM_CXXFLAGS@

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
#noinst_LIBRARIES = libimageClasses.a
lib_LTLIBRARIES = libimageClasses.la

# where to install the headers on the system
libimageClasses_ladir = $(includedir)/pktools/imageclasses

## Instruct libtool to include ABI version information in the generated shared
## library file (.so).  The library ABI version is defined in configure.ac, so
## that all version information is kept in one place.
libimageClasses_la_LDFLAGS = -version-info $(PKTOOLS_SO_VERSION) $(AM_LDFLAGS)

# the list of header files that belong to the library (to be installed later)
libimageClasses_la_HEADERS = ImgReaderGdal.h  ImgReaderOgr.h  ImgWriterGdal.h  ImgWriterOgr.h

# the sources to add to the library and to add to the source distribution
libimageClasses_la_SOURCES = $(libimageClasses_la_HEADERS) ImgReaderGdal.cc  ImgReaderOgr.cc  ImgWriterGdal.cc  ImgWriterOgr.cc
###############################################################################
