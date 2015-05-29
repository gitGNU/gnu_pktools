# to include m4 macros into aclocal.m4
ACLOCAL_AMFLAGS = -I m4

# the subdirectories of the project to go into
SUBDIRS = src/base \
	src/imageclasses \
	src/algorithms \
	src/fileclasses \
	$(LASCLASSES_OPT) \
	$(DOX_OPT) \
	src/apps

## The generated configuration header is installed in its own subdirectory of
## $(libdir).  The reason for this is that the configuration information put
## into this header file describes the target platform the installed library
## has been built for.  Thus the file must not be installed into a location
## intended for architecture-independent files, as defined by the Filesystem
## Hierarchy Standard (FHS).
## The nodist_ prefix instructs Automake to not generate rules for including
## the listed files in the distribution on 'make dist'.  Files that are listed
## in _HEADERS variables are normally included in the distribution, but the
## configuration header file is generated at configure time and should not be
## shipped with the source tarball.
pktools_libincludedir = $(libdir)/pktools/include
nodist_pktools_libinclude_HEADERS = config.h


## Install the generated pkg-config file (.pc) into the expected location for
## architecture-dependent package configuration information.  Occasionally,
## pkg-config files are also used for architecture-independent data packages,
## in which case the correct install location would be $(datadir)/pkgconfig.
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pktools.pc


########################################################################
#       Additional installation for QGIS processing plugin files       #
########################################################################
qgisdir = $(HOME)/.qgis2/python/plugins/$(PACKAGE)
qgis_DATA = \
	qgis/pktools/logo.png \
	qgis/pktools/__init__.py \
	qgis/pktools/metadata.txt \
	qgis/pktools/ProcessingPktoolsPlugin.py \
	qgis/pktools/pktoolsAlgorithmProvider.py \
	qgis/pktools/pktoolsAlgorithm.py \
	qgis/pktools/pktoolsUtils.py \
	qgis/pktools/pkcomposite.py \
	qgis/pktools/pkcrop.py \
	qgis/pktools/pkreclass.py \
	qgis/pktools/pkextract_grid.py \
	qgis/pktools/pkextract.py \
	qgis/pktools/pkextract_random.py \
	qgis/pktools/pkfilterdem.py \
	qgis/pktools/pkfilter_spatial.py \
	qgis/pktools/pkfilter_spectral.py \
	qgis/pktools/pkgetmask.py \
	qgis/pktools/pklas2img.py \
	qgis/pktools/pksetmask.py \
	qgis/pktools/pksvm.py

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = doc $(qgis_DATA)
#bin_SCRIPTS = pktools-config
