include(CTest)

set(UNIT_TESTS
    test_affine
    test_affine_model
    test_function_operations
    test_box_function
    test_function_dispatch
    test_function_interface
    test_function_models
    test_symbolic_function
    test_function_names
    test_function_patch_operations
    test_projection
    test_domain_scaling
    test_taylor_series
    test_function
    test_formula
    test_measurable_function
    test_multifunction
    test_procedure
    test_scaled_function_patch
    test_taylor_function
    test_taylor_model
)
foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} PRIVATE ariadne-function function)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "function")
    if(WIN32)
        set_tests_properties(${TEST} PROPERTIES
            ENVIRONMENT_MODIFICATION
                "PATH=path_list_prepend:$<TARGET_FILE_DIR:ariadne-function>"
        )
    endif()
endforeach()
add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS})
