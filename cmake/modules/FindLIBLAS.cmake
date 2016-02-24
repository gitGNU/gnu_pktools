
find_path(LIBLAS_INCLUDE_DIR liblas/liblas.hpp
    HINTS
        ENV LIBLAS_DIR
        ENV LIBLAS_ROOT
    PATH_SUFFIXES
        include
)

find_library(LIBLAS_LIBRARY las
    HINTS
        ENV LIBLAS_DIR
        ENV LIBLAS_ROOT
    PATH_SUFFIXES
        lib
        liblas/
)

find_library(LIBLAS_C_LIBRARY las_c
    HINTS
        ENV LIBLAS_DIR
        ENV LIBLAS_ROOT
    PATH_SUFFIXES
        lib
        liblas/
)

set(LIBLAS_INCLUDE_DIRS "${LIBLAS_INCLUDE_DIR}")
set(LIBLAS_LIBRARIES "${LIBLAS_LIBRARY}" "${LIBLAS_C_LIBRARY}")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LIBLAS DEFAULT_MSG
    LIBLAS_LIBRARY LIBLAS_INCLUDE_DIR LIBLAS_C_LIBRARY)

mark_as_advanced(LIBLAS_INCLUDE_DIR LIBLAS_LIBRARY LIBLAS_C_LIBRARY)
