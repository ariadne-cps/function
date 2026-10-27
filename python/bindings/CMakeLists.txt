include_directories(${Python3_INCLUDE_DIRS})

get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-function-submodule-obj OBJECT
    function_submodule.cpp
)
target_include_directories(pyariadne-function-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    $<TARGET_PROPERTY:pyariadne-geometry,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-function-submodule-obj PRIVATE
    pyariadne-python-common
    pybind11::module
    function
)
set_target_properties(pyariadne-function-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
    target_compile_options(pyariadne-function-submodule-obj PRIVATE
        /bigobj
        /wd4702
    )
endif()

add_library(pyariadne-function INTERFACE)
target_sources(pyariadne-function INTERFACE
    $<TARGET_OBJECTS:pyariadne-function-submodule-obj>
)
target_include_directories(pyariadne-function INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-geometry,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-function INTERFACE
    pyariadne-geometry
)

ariadne_add_python_module(
    NAME function
    MODULE_SOURCE function_module.cpp
)
