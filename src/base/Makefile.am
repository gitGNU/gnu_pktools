###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
#do not want to install pktestoption
noinst_PROGRAMS = pktestOption  

# additional include pathes necessary to compile the C++ library
AM_LDFLAGS = $(GDAL_LDFLAGS) @AM_LDFLAGS@
AM_CXXFLAGS = -I$(top_srcdir)/src $(GDAL_CFLAGS) @AM_CXXFLAGS@
#AM_CXXFLAGS = -I$(top_srcdir)/src @AM_CXXFLAGS@

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

#noinst_LIBRARIES = libbase.a
lib_LTLIBRARIES = libbase.la

# where to install the headers on the system
libbase_ladir = $(includedir)/pktools/base

## Instruct libtool to include ABI version information in the generated shared
## library file (.so).  The library ABI version is defined in configure.ac, so
## that all version information is kept in one place.
libbase_la_LDFLAGS = -version-info $(PKTOOLS_SO_VERSION) $(AM_LDFLAGS)

# the list of header files that belong to the library (to be installed later)
libbase_la_HEADERS = Vector2d.h IndexValue.h Optionpk.h PosValue.h

# the sources to add to the library and to add to the source distribution
###############################################################################
libbase_la_SOURCES = $(libbase_la_HEADERS)

# list of sources for the binaries
pktestOption_SOURCES = pktestOption.cc

###############################################################################
# HEADERS USED FOR BUILDING BUT NOT TO BE INSTALLED
###############################################################################

#noinst_HEADERS = Optionpk.h PosValue.h Vector2d.h
###############################################################################
