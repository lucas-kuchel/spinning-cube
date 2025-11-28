cmake_minimum_required(VERSION 3.21)

project(engine LANGUAGES C)

file(GLOB_RECURSE SOURCES "src/*.c")

add_executable(cube ${SOURCES})

target_include_directories(cube PRIVATE
    "include"
)

target_link_libraries(cube PRIVATE m)
