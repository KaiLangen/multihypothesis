cmake_minimum_required(VERSION 2.8)

project( obme )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
configure_file(config.h.in ${CMAKE_BINARY_DIR}/generated/config.h)
include_directories(${CMAKE_BINARY_DIR}/generated/)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -da")
find_package(OpenCV REQUIRED)
#find_package( OpenCV REQUIRED )

#execute_process (
#  COMMAND bash -c "cp ${CMAKE_CURRENT_SOURCE_DIR}/bowing.cfg ${CMAKE_BINARY_DIR}" 
#)

execute_process (
  COMMAND bash -c "tar -xzvf ${CMAKE_CURRENT_SOURCE_DIR}/ldpca.tar.gz \
                   --directory ${CMAKE_BINARY_DIR}"
)

execute_process (
  COMMAND bash -c "tar -xzvf ${CMAKE_CURRENT_SOURCE_DIR}/jm.tar.gz \
                   --directory ${CMAKE_BINARY_DIR}"
)

add_library(utils bitstream.cpp
                  cavlc.cpp
                  codec.cpp
                  frameBuffer.cpp
                  ldpca.cpp
                  transform.cpp
                  calculations.cpp
)

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
                       sideInformation_MCI.cpp
                       sideInformation_chromaME.cpp
)
target_link_libraries(decoder utils m)

# Tester
add_executable(tester  encoder.cpp
                       cavlcEnc.cpp
                       cavlcDec.cpp
                       ldpcaEnc.cpp
                       chroma_test.cpp
)

target_link_libraries(tester utils m ${OpenCV_LIBS})
