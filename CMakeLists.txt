cmake_minimum_required(VERSION 3.17)
project(ssdg)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(lib)
add_subdirectory(catch)
add_executable(test catch_main.cpp test.cpp)


target_link_libraries(test sim_annealing_lib catch)


