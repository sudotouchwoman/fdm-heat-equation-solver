cmake_minimum_required(VERSION 3.16)
project(finite-diff-method)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Werror -Wpedantic)

set(PROJECT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(project)
