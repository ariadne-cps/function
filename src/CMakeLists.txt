add_library(FUNCTION_SRC OBJECT
    box_function.cpp
    affine.cpp
    formula.cpp
    procedure.cpp
    function.cpp
    measurable_function.cpp
    multifunction.cpp
    function_patch.cpp
    function_model.cpp
    affine_model.cpp
    taylor_model.cpp
    taylor_function.cpp
    taylor_multifunction.cpp
    symbolic_function.cpp
)
target_include_directories(FUNCTION_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/function"
)
