cmake_minimum_required(VERSION 3.15)

project(Function VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT TARGET ariadne-configuration)
    find_package(Git REQUIRED)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}"
                    submodule update --init -- ss/cfg
                    RESULT_VARIABLE _result)
    if(_result)
        message(FATAL_ERROR "Cannot initialize configuration submodule")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/ss/cfg/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        CHAINS
            "configuration"
            "geometry|configuration"
    )
endif()

setup_project()

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(ss)

add_library(function INTERFACE)
target_include_directories(function INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-function>
)
target_link_libraries(function INTERFACE geometry)

add_subdirectory(src)
target_link_libraries(FUNCTION_SRC PRIVATE function)

ariadne_register_public_headers(
    TARGET function
    DIRECTORIES "${PROJECT_SOURCE_DIR}/include/function"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-function.hpp"
)

setup_project_library(
    TARGET ariadne-function
    OBJECTS FUNCTION_SRC
    LINK_LIBRARIES function
)

setup_standalone_project_tests(
    TARGET ariadne-function
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    ariadne_install_dependency_bundle(
        TARGET function
        DESTINATION include/ariadne-function
    )
endif()
