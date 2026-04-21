cmake_minimum_required(VERSION 3.10)

project(PlotExample)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python3 COMPONENTS Development REQUIRED)

# Get NumPy include directory
execute_process(
    COMMAND python3 -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${Python3_INCLUDE_DIRS}
    ${NUMPY_INCLUDE_DIR}
)

# Source files
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/source/*.cc)

add_executable(plot_app ${SOURCES})

target_link_libraries(plot_app PRIVATE Python3::Python)