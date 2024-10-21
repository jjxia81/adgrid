cmake_minimum_required(VERSION 3.14)

project(adaptive_mesh_refinement LANGUAGES CXX C)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(cli11)
include(Eigen3)
include(nlohmann-json)

# project library
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.c* ${CMAKE_CURRENT_LIST_DIR}/src/3rd/*.c* 
                    ${CMAKE_CURRENT_LIST_DIR}/src/3rd/implicit_functions/*.c* ${CMAKE_CURRENT_LIST_DIR}/src/*.h 
                    ${CMAKE_CURRENT_LIST_DIR}/src/3rd/mtet/*.c* ${CMAKE_CURRENT_LIST_DIR}/src/3rd/mshio/*.c*
                    ${CMAKE_CURRENT_LIST_DIR}/src/3rd/nanothread/*.c*)
add_library(adaptive_mesh_refinement STATIC ${SRC_FILES})
target_include_directories(adaptive_mesh_refinement PUBLIC
     ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(adaptive_mesh_refinement
    PUBLIC 
    Eigen3::Eigen
    nlohmann_json::nlohmann_json)

add_executable(gridgen "app/gridgen")

target_compile_features(gridgen PRIVATE cxx_std_20)
target_include_directories(gridgen PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(gridgen PRIVATE
    adaptive_mesh_refinement
    CLI11::CLI11
    )

option(GRID_GEN_TESTS "Build adaptive grid gen unit tests" ON)
if (GRID_GEN_TESTS)
 include(CTest)
 enable_testing()
 include(Catch2)

 file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
 add_executable(grid_gen_tests ${TEST_FILES})
 target_link_libraries(grid_gen_tests
    PRIVATE 
    adaptive_mesh_refinement
    Catch2::Catch2)
 target_compile_features(grid_gen_tests PRIVATE cxx_std_20)
 target_compile_definitions(grid_gen_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
target_compile_definitions(
    grid_gen_tests
    PRIVATE                   
    -DTEST_FILE="${CMAKE_CURRENT_LIST_DIR}/data"
)
 catch_discover_tests(grid_gen_tests)
endif()
