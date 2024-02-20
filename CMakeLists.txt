cmake_minimum_required (VERSION 3.5)
project(MyRayTracer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# file(GLOB SOURCES ./src/*.h ./src/*.cpp)
file(GLOB HEADER_FILES ./src/*.h)
file(GLOB SOURCE_FILES ./src/*.cpp)

add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES})

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
