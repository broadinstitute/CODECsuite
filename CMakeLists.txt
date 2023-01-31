cmake_minimum_required(VERSION 3.18)
set( CMAKE_CXX_STANDARD 14 )
project (codecsuite LANGUAGES CXX)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
#for gprof
#set (CMAKE_CXX_FLAGS "${CAMKE_CXX_FLAGS} -g -pg")
#for valgrind
#set (CMAKE_CXX_FLAGS "${CAMKE_CXX_FLAGS} -g -O1")
set( CMAKE_EXE_LINKER_FLAGS " ${CMAKE_EXE_LINKER_FLAGS} -pthread")
set( CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -Wno-unused-variable " )


set (CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -DSEQAN_ENABLE_TESTING=0 -DSEQAN_HAS_ZLIB=1")
find_package(OpenMP REQUIRED)

include_directories(include)
include_directories(bbcpputil/third_party/seqan/include/)

set( CMAKE_EXE_LINKER_FLAGS " ${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -L. ")
set( CPP_LIBS ${CPP_LIBS} seqlib bwa fml hts bz2)
set( CPP_LIBS ${CPP_LIBS} dl )
set( CPP_LIBS ${CPP_LIBS} z lzma)
set( CPP_LIBS ${CPP_LIBS} OpenMP::OpenMP_CXX)

include_directories(bbcpputil/include)

set(SEQLIB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/bbcpputil/third_party/SeqLib)
find_library(SEQLIB_LIBRARY
        NAMES seqlib
        PATHS ${SEQLIB_ROOT}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH
        )
include(ExternalProject)
if(NOT SEQLIB_LIBRARY)
    message("SeqLib not found")
    ExternalProject_Add(SeqLib
            SOURCE_DIR "${SEQLIB_ROOT}"
            INSTALL_DIR "${SEQLIB_ROOT}"
            BUILD_IN_SOURCE 1
            CONFIGURE_COMMAND chmod +x ./configure COMMAND ./configure prefix=${SEQLIB_ROOT}
            BUILD_COMMAND make CXXFLAGS='-std=c++11'
            INSTALL_COMMAND make install
            )
else()
    message("SeqLib found")
    add_custom_target(SeqLib)
endif()

include_directories(${SEQLIB_ROOT})
include_directories(${SEQLIB_ROOT}/htslib)
link_directories(${SEQLIB_ROOT}/lib/)

## not currently used. But reserve for denovo consensus
if (NOT TARGET spoa)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bbcpputil/third_party/spoa EXCLUDE_FROM_ALL)
endif ()

set(deps bbcpputil/cpp/AlignmentConsensus.cpp bbcpputil/cpp/BamRecordExt.cpp)
set(subcommand demux.cpp trimadapter.cpp consensus.cpp accuracy.cpp)
if (CMAKE_BUILD_TYPE MATCHES Debug)
    message(" ")
    message("CMAKE IN DEBUG MODE")
    message(" ")
    link_directories($ENV{ELFUTILS_ROOT}/lib/)
    set(LIBDW_LIBRARY "$ENV{ELFUTILS_ROOT}/lib/libdw.a")
    set(LIBDW_INCLUDE_DIR "$ENV{ELFUTILS_ROOT}/include/")
    set( CPP_LIBS ${CPP_LIBS} dw elf) # backward_cpp
    add_subdirectory(bbcpputil/third_party/backward-cpp)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    message(" ")
    message("CMAKE IN RELEASE MODE")
    message(" ")
endif ()

add_executable(codec codec.cpp ${subcommand} ${deps})
add_dependencies(codec SeqLib)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    #add_executable(codec codec.cpp ${subcommand} ${deps} ${BACKWARD_ENABLE})
    add_backward(codec)
endif()
target_link_libraries(codec ${CPP_LIBS} spoa)