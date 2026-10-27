find_library(ARIADNE_FUNCTION_LIBRARY NAMES ariadne-function)
find_path(ARIADNE_FUNCTION_INCLUDE_DIR ariadne-function.hpp PATH_SUFFIXES ariadne-function)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneFunction DEFAULT_MSG
  ARIADNE_FUNCTION_LIBRARY
  ARIADNE_FUNCTION_INCLUDE_DIR
)

if(AriadneFunction_FOUND)
  get_filename_component(ARIADNE_FUNCTION_INCLUDE_PARENT_DIR ${ARIADNE_FUNCTION_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_FUNCTION_INCLUDE_DIRS
    ${ARIADNE_FUNCTION_INCLUDE_PARENT_DIR}
    ${ARIADNE_FUNCTION_INCLUDE_DIR}
  )
  if(NOT TARGET AriadneFunction::ariadne-function)
    add_library(AriadneFunction::ariadne-function UNKNOWN IMPORTED)
    set_target_properties(AriadneFunction::ariadne-function PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_FUNCTION_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_FUNCTION_INCLUDE_DIRS}"
    )
  endif()
endif()
