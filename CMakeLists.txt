cmake_minimum_required(VERSION 3.21)

project(PROJ-MO824
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Set a default build type if none was specified
set(default_build_type "Release")

# Headers directories
include_directories("include")
include_directories("lemon/lemon-1.3.1/include")
include_directories("brkga/src")
include_directories($ENV{GUROBI_HOME}/include)

# Utils used in all targets
add_library(graph_utils src/mygraphlib.cpp src/geompack.cpp src/myutils.cpp
        src/mycolor.cpp src/pickup_delivery_utils.cpp)

# Lemon lib directory
link_directories(${CMAKE_SOURCE_DIR}/lemon/lemon-1.3.1/lib)

# Gurobi lib directory
link_directories($ENV{GUROBI_HOME}/lib)

# Compiler flags
set(CMAKE_CXX_FLAGS "-m64 -Ofast -Wall -Wunused-local-typedefs -lemon -fopenmp")

# Target executable and linking
add_executable(pickup_delivery_branch_and_bound src/pickup_delivery_branch_and_bound.cpp)
target_link_libraries(pickup_delivery_branch_and_bound m)
target_link_libraries(pickup_delivery_branch_and_bound graph_utils)

add_executable(pickup_delivery_heuristic src/pickup_delivery_heuristic.cpp)
target_link_libraries(pickup_delivery_heuristic m)
target_link_libraries(pickup_delivery_heuristic graph_utils)

add_executable(pickup_delivery_brkga src/pickup_delivery_brkga.cpp src/PickupDeliveryDecoder.cpp)
target_link_libraries(pickup_delivery_brkga m)
target_link_libraries(pickup_delivery_brkga pthread)
target_link_libraries(pickup_delivery_brkga graph_utils)

add_executable(pickup_delivery_ilp src/pickup_delivery_ilp.cpp)
target_link_libraries(pickup_delivery_ilp m)
target_link_libraries(pickup_delivery_ilp pthread)
target_link_libraries(pickup_delivery_ilp graph_utils)
target_link_libraries(pickup_delivery_ilp gurobi_c++ gurobi95)
