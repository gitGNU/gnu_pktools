AC_INIT([pktools], [2.6.2], [kempenep@gmail.com])
#AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_INIT_AUTOMAKE([-Wall -Wno-extra-portability foreign])
#AM_INIT_AUTOMAKE([subdir-objects]) #not working due to bug in autoconf, see Debian list: Bug #752993)
AC_CONFIG_MACRO_DIR([m4])

AX_LIB_GDAL()
GDAL_FOUND="no"
if test ! -z "$GDAL_CFLAGS" -a ! -z "$GDAL_LDFLAGS"; then
GDAL_FOUND="yes"
fi
dnl AC_MSG_WARN([$GDAL_VERSION])
dnl AC_MSG_WARN([$GDAL_CFLAGS])
dnl AC_MSG_WARN([$GDAL_LDFLAGS])
dnl AC_MSG_WARN([$GDAL_DEP_LDFLAGS])
dnl AC_MSG_WARN([$GDAL_OGR_ENABLED])

# check for C++ compiler and the library compiler
AC_PROG_CXX
AC_PROG_FC
LT_INIT

# check if the source folder is correct
AC_CONFIG_SRCDIR([src/apps/pkinfo.cc])
AX_LIB_GDAL([1.10.0]) dnl uncomment if gdal version 1.10.0 is required

AC_CHECK_HEADERS([gdal.h])

dnl Default-disabled features
AC_ARG_ENABLE([fann],
AS_HELP_STRING([--enable-fann], [Enable feature fann]))

AS_IF([test "x$enable_fann" = "xyes"], [
PKG_CHECK_MODULES(FANN, fann >= 2.1.0, [HAVE_FANN=1], [HAVE_FANN=0])
AM_CONDITIONAL([USE_FANN], [test "$HAVE_FANN" -eq 1])
dnl AC_MSG_WARN("USE_FANN=true")
], [
AM_CONDITIONAL(USE_FANN, false)
dnl AC_MSG_WARN("USE_FANN=false")
])

dnl AC_ARG_ENABLE([svm],
dnl     AS_HELP_STRING([--enable-svm], [Enable feature svm]))

dnl AS_IF([test "x$enable_svm" = "xyes"], [
dnl 	AM_CONDITIONAL(USE_SVM, true)
dnl 	AC_MSG_WARN("USE_SVM=true")
dnl 	])

AC_ARG_ENABLE([las],
AS_HELP_STRING([--enable-las], [Enable feature las]))

AS_IF([test "x$enable_las" = "xyes"], [
AM_CONDITIONAL(USE_LAS, true)
dnl AC_MSG_WARN("USE_LAS=true")
LASCLASSES_OPT="src/lasclasses"
AC_SUBST([LASCLASSES_OPT])
], [
AM_CONDITIONAL(USE_LAS, false)
dnl AC_MSG_WARN("USE_LAS=false")
])

AC_ARG_ENABLE([nlopt],
AS_HELP_STRING([--enable-nlopt], [Enable feature nlopt]))

AS_IF([test "x$enable_nlopt" = "xyes"], [
PKG_CHECK_MODULES(NLOPT, nlopt >= 2.1.0, [HAVE_NLOPT=1], [HAVE_NLOPT=0])
AM_CONDITIONAL([USE_NLOPT], [test "$HAVE_NLOPT" -eq 1])
dnl AC_MSG_WARN("USE_NLOPT=true")
], [
AM_CONDITIONAL(USE_NLOPT, false)
dnl AC_MSG_WARN("USE_NLOPT=false")
])


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

# Define these substitions here to keep all version information in one place.
# For information on how to properly maintain the library version information,
# refer to the libtool manual, section "Updating library version information":
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
AC_SUBST([PKTOOLS_SO_VERSION], [1:3:0])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_HEADERS([config.h])
#frankie says: "You should not distribute the config.h file, but the config.h.in file only".
#AC_CONFIG_HEADERS([config.h])

AM_COND_IF([USE_LAS],
	[AC_CONFIG_FILES([
	src/lasclasses/Makefile])])
AC_CONFIG_FILES([
	Makefile
	src/base/Makefile
	src/algorithms/Makefile
	src/imageclasses/Makefile
	src/fileclasses/Makefile
	src/apps/Makefile
	pktools.pc:pktools.pc.in
	])
# generate the final Makefile etc.
AC_OUTPUT
