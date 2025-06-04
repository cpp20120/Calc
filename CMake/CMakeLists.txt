cmake_minimum_required(VERSION 3.28)
project(Calc)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED  ON)
#set(CMAKE_EXPORT_COMPILE_COMMAND ON)

add_executable(${CMAKE_PROJECT_NAME} src/main.cpp)
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE include)
