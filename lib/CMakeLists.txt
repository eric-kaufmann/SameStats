cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 14)


add_library(sim_annealing_lib STATIC include/sim_annealing.h src/sim_annealing.cpp include/utils.h src/utils.cpp include/line_shapes.h src/line_shapes.cpp)
target_include_directories(sim_annealing_lib PUBLIC include PRIVATE src)

find_package(OpenMP REQUIRED)
target_link_libraries(sim_annealing_lib PUBLIC OpenMP::OpenMP_CXX)
