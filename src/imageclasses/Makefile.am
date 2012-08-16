AM_LDFLAGS = $(GDAL_LDFLAGS) @AM_LDFLAGS@
AM_CXXFLAGS = -I$(top_srcdir)/src $(GDAL_CFLAGS) @AM_CXXFLAGS@

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
noinst_LIBRARIES = libimageClasses.a

# where to install the headers on the system
libimageClasses_adir = $(includedir)/imageclasses

# the list of header files that belong to the library (to be installed later)
libimageClasses_a_HEADERS = ImgReaderGdal.h  ImgReaderOgr.h  ImgWriterGdal.h  ImgWriterOgr.h

# the sources to add to the library and to add to the source distribution
libimageClasses_a_SOURCES = $(libimageClasses_a_HEADERS) ImgReaderGdal.cc  ImgReaderOgr.cc  ImgWriterGdal.cc  ImgWriterOgr.cc
###############################################################################
