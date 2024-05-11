cmake_minimum_required(VERSION 3.10)

file(GLOB SRCS *.cpp)

add_executable(lab3 ${SRCS})

target_include_directories(lab3
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../eigen-3.4.0/Eigen
)
