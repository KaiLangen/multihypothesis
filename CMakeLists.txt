cmake_minimum_required(VERSION 2.8)

project( obme )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -da")
find_package( OpenCV REQUIRED )
add_executable(obme decoder.cpp
                    lowpass.cpp
                    obmc.cpp
                    obme.cpp
                    obme.h
                    resampling.cpp
                    types.h
)

target_link_libraries(obme ${OpenCV_LIBS})


