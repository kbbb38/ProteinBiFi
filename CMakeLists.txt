cmake_minimum_required(VERSION 3.10)

project(ProteinBiFi LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_ASAN "Enable AddressSanitizer" OFF)

if(ENABLE_ASAN)
    set(SANITIZER_FLAGS "-fsanitize=address -fno-omit-frame-pointer -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZER_FLAGS}")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Ofast -pthread")

add_executable(ProteinBiFi
    ProteinBiFi.cpp
    src/SpectrumBitSet.h
    src/SpectrumBitSet.cpp
    src/ExperimentalBitSets.h
    src/ExperimentalBitSets.cpp
    src/Constants.h
)