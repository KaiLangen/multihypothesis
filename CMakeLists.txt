cmake_minimum_required(VERSION 2.8)

project( obme )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
configure_file(config.h.in ${CMAKE_BINARY_DIR}/generated/config.h)
include_directories(${CMAKE_BINARY_DIR}/generated/)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -da")

execute_process (
  COMMAND bash -c "tar -xzvf ${CMAKE_CURRENT_SOURCE_DIR}/ldpca.tar.gz \
                   --directory ${CMAKE_BINARY_DIR}"
)

execute_process (
  COMMAND bash -c "tar -xzvf ${CMAKE_CURRENT_SOURCE_DIR}/jm.tar.gz \
                   --directory ${CMAKE_BINARY_DIR}"
)

execute_process (
  COMMAND bash -c "if [ `uname -m` == x86_64 ]
  then
    mv ${CMAKE_BINARY_DIR}/jm/lencod_64.exe ${CMAKE_BINARY_DIR}/jm/lencod.exe
  elif [ `uname -m` == x86_64 ]
  then
    mv ${CMAKE_BINARY_DIR}/jm/lencod_32.exe ${CMAKE_BINARY_DIR}/jm/lencod.exe
  fi"
)

add_library(utils bitstream.cpp
                  cavlc.cpp
                  codec.cpp
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
                       sideInformation_MCI.cpp
                       sideInformation_chromaME.cpp
                       dec_main.cpp
)
target_link_libraries(decoder utils m)

# Tester
add_executable(tester  encoder.cpp
                       cavlcEnc.cpp
                       cavlcDec.cpp
                       ldpcaEnc.cpp
                       chroma_test.cpp
)

target_link_libraries(tester utils m)
