AC_INIT([pktools], [2.1], [kempenep@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
AX_LIB_GDAL()
GDAL_FOUND="no"
if test ! -z "$GDAL_CFLAGS" -a ! -z "$GDAL_LDFLAGS"; then
    GDAL_FOUND="yes"
fi
AC_MSG_WARN([$GDAL_VERSION])
AC_MSG_WARN([$GDAL_CFLAGS])
AC_MSG_WARN([$GDAL_LDFLAGS])
AC_MSG_WARN([$GDAL_DEP_LDFLAGS])
AC_MSG_WARN([$GDAL_OGR_ENABLED])

# check for C++ compiler and the library compiler
AC_PROG_CXX
AC_PROG_RANLIB

# check if the source folder is correct
AC_CONFIG_SRCDIR([src/apps/pkinfo.cc])
AX_LIB_GDAL([1.8.0]) dnl uncomment if gdal version 1.8.0 is required

AC_CHECK_HEADERS([gdal.h])

# Checks for packages
#PKG_CHECK_MODULES(FANN, fann >= 2.1.0, [HAVE_FANN=1], [HAVE_FANN=0])
#AM_CONDITIONAL([USE_FANN], [test "$HAVE_FANN" -eq 1])

AC_ARG_WITH([fann],
	[AS_HELP_STRING([--with-fann], [use fann library] @@)],
	[],
	[with_fann=no])
AS_IF([test "x$with_fann" != xno],
	[
	PKG_CHECK_MODULES(FANN, fann >= 2.1.0, [HAVE_FANN=1], [HAVE_FANN=0])
	AM_CONDITIONAL([USE_FANN], [test "$HAVE_FANN" -eq 1])
	], [
	AM_CONDITIONAL(USE_FANN, false)
	])

AC_ARG_WITH([las],
	[AS_HELP_STRING([--with-las], [use las library] @@)],
	[],
	[with_las=no])
AS_IF([test "x$with_las" != xno],
	[
	AM_CONDITIONAL(USE_LAS, true)
	], [
	AM_CONDITIONAL(USE_LAS, false)
	])

if $USE_LAS; then
   AC_MSG_WARN("USE_LAS=true")
   FILECLASSES_OPT="src/fileclasses"
   AC_SUBST([FILECLASSES_OPT])
fi

PKG_CHECK_MODULES(GSL, gsl >= 1.14, [HAVE_GSL=1], [HAVE_GSL=0])
AM_CONDITIONAL([USE_GSL], [test "$HAVE_GSL" -eq 1])

# use the C++ compiler for the following checks
AC_LANG([C++])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
src/base/Makefile
src/algorithms/Makefile
src/imageclasses/Makefile
src/fileclasses/Makefile
src/apps/Makefile
])

# generate the final Makefile etc.
AC_OUTPUT
