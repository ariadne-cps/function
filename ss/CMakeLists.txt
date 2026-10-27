if(NOT TARGET geometry)
    ensure_submodule_initialized(geometry)
    add_subdirectory(geometry)
endif()
