cmake_minimum_required (VERSION 3.5)
project(MyRayTracer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES ./src/*.h ./src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})


