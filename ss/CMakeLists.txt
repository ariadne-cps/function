if(NOT TARGET geometry)
    ariadne_add_submodule(geometry)
endif()

# Testing remains a test-only dependency.
ariadne_setup_testing()
