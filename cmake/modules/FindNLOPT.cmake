find_path(NLOPT_INCLUDE_DIR nlopt.h
    HINTS
        ENV NLOPT_DIR
        ENV NLOPT_ROOT
    PATH_SUFFIXES
        include
        include/nlopt
)

find_library(NLOPT_LIBRARY libnlopt-0
    HINTS
        ENV NLOPT_DIR
        ENV NLOPT_ROOT
    PATH_SUFFIXES
        lib
        lib/nlopt
)

set(NLOPT_INCLUDE_DIRS "${NLOPT_INCLUDE_DIR}")
set(NLOPT_LIBRARIES "${NLOPT_LIBRARY}")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(NLOPT DEFAULT_MSG
    NLOPT_LIBRARY NLOPT_INCLUDE_DIR)

mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARY)