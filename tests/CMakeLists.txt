include(CTest)
include(Catch)

foreach (CDIM 2)
    set(TEST_EXEC_NAME "tests_cadmus_${CDIM}")
    add_executable(${TEST_EXEC_NAME} test.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/coh_block.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/cube.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/cubical_filtration.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/filtration.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/lower_star.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/matrix.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/matrix_with_values.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/messages.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/simplex.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/types.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../src/values_only_matrix.cpp
            )

    target_compile_definitions(${TEST_EXEC_NAME} PRIVATE DIY_USE_SPDLOG)
    target_compile_definitions(${TEST_EXEC_NAME} PRIVATE CADMUS_DIM=2)
    target_compile_options(${TEST_EXEC_NAME} PRIVATE -Wno-deprecated)
    target_link_libraries(${TEST_EXEC_NAME} PRIVATE ${libraries})
    target_link_libraries(${TEST_EXEC_NAME} PRIVATE Catch2::Catch2WithMain)

    catch_discover_tests(${TEST_EXEC_NAME})
endforeach()