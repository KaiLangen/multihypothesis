cmake_minimum_required(VERSION 2.8)

project( obme )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
configure_file(config.h.in ${CMAKE_BINARY_DIR}/generated/config.h)
include_directories(${CMAKE_BINARY_DIR}/generated/)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -da")
find_package( OpenCV REQUIRED )

execute_process (
  COMMAND bash -c "date"
  COMMAND bash -c "tar -xzvf ${CMAKE_CURRENT_SOURCE_DIR}/ldpca.tar.gz --directory ${CMAKE_BINARY_DIR}"
  OUTPUT_VARIABLE outVar
)

message (
  STATUS "${outVar}"
)

add_library(utils bitstream.cpp
                  cavlc.cpp
                  codec.cpp
                  ldpca.cpp
                  transform.cpp
)

# test bin
#add_executable(obme decoder.cpp
#                    lowpass.cpp
#                    obmc.cpp
#                    obme.cpp
#                    obme.h
#                    resampling.cpp
#                    types.h
#)
#target_link_libraries(obme ${OpenCV_LIBS})

# encoder
add_executable(encoder encoder.cpp
                       cavlcEnc.cpp
                       ldpcaEnc.cpp
                       enc_main.cpp
)
target_link_libraries(encoder utils)


# decoder
add_executable(decoder decoder.cpp
                       cavlcDec.cpp
                       corrModel.cpp
                       ldpcaDec.cpp
                       dec_main.cpp
                       sideInformation.cpp
)
target_link_libraries(decoder utils ${OpenCV_LIBS})

