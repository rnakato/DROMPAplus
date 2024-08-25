cmake_minimum_required(VERSION 2.8)
project(DROMPAPLUS CXX)

set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_CXX_FLAGS "-std=c++14 -O2 -W -Wall")
if(ENABLE_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG") # -fconcepts
endif()

include_directories("/usr/local/include")

add_subdirectory(src)
add_subdirectory(submodules/SSP/)
add_subdirectory(test)
