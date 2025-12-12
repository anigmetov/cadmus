foreach (CDIM 1 2 3)
    set(EXECNAME "cadmus_${CDIM}")
    add_executable(${EXECNAME}
            main.cpp
            simplex.cpp
            cube.cpp
            filtration.cpp
            cubical_filtration.cpp
            matrix.cpp
            matrix_with_values.cpp
            values_only_matrix.cpp
            lower_star.cpp
            messages.cpp
            coh_block.cpp
            types.cpp)

    target_compile_definitions(${EXECNAME} PRIVATE DIY_USE_SPDLOG)
    target_compile_definitions(${EXECNAME} PRIVATE CADMUS_DIM=${CDIM})
    target_compile_options(${EXECNAME} PRIVATE -Wno-deprecated)
    target_link_libraries(${EXECNAME} PRIVATE ${libraries})
    if (cadmus_backward)
        target_link_libraries(${EXECNAME} PUBLIC Backward::Object)
    endif()
endforeach ()

