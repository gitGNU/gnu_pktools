###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
#do not want to install pktestoption
noinst_PROGRAMS = pktestOption  

# additional include pathes necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/src @AM_CXXFLAGS@

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
noinst_LIBRARIES = libbase.a

# where to install the headers on the system
libbase_adir = $(includedir)/base

# the list of header files that belong to the library (to be installed later)
libbase_a_HEADERS = $(top_srcdir)/config.h PointData.h Vector2d.h IndexValue.h Optionpk.h PosValue.h

# the sources to add to the library and to add to the source distribution
libbase_a_SOURCES = $(libbase_a_HEADERS) PointData.cc
###############################################################################

# list of sources for the binaries
pktestOption_SOURCES = pktestOption.cc

###############################################################################
# HEADERS USED FOR BUILDING BUT NOT TO BE INSTALLED
###############################################################################

noinst_HEADERS = Optionpk.h PosValue.h Vector2d.h PointData.h
###############################################################################
