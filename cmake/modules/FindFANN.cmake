
find_path(FANN_INCLUDE_DIR fann.h
    HINTS
        ENV FANN_DIR
        ENV FANN_ROOT
    PATH_SUFFIXES
        include
        include/fann
)

find_library(FANN_LIBRARY floatfann
    HINTS
        ENV FANN_DIR
        ENV FANN_ROOT
    PATH_SUFFIXES
        lib
        lib/fann
)

set(FANN_INCLUDE_DIRS "${FANN_INCLUDE_DIR}")
set(FANN_LIBRARIES "${FANN_LIBRARY}")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FANN DEFAULT_MSG
    FANN_LIBRARY FANN_INCLUDE_DIR)

mark_as_advanced(FANN_INCLUDE_DIR FANN_LIBRARY)
