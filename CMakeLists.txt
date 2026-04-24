cmake_minimum_required(VERSION 3.10)
project(simDER)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-fopenmp)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

set(MKLROOT /opt/intel/oneapi/mkl/latest)

include_directories(
    ${EIGEN3_INCLUDE_DIR}
    ${MKLROOT}/include
)

link_directories(
    ${MKLROOT}/lib/intel64
)

file(GLOB SOURCES "*.cpp")

add_executable(simDER ${SOURCES})

target_link_libraries(simDER
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    -Wl,--start-group
    mkl_intel_lp64
    mkl_intel_thread
    mkl_core
    -Wl,--end-group
    iomp5
    pthread
    m
    dl   # 很多时候漏掉 dl 库会导致 MKL 加载失败
)