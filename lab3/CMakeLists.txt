cmake_minimum_required(VERSION 2.8)

project(Lab3)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(max_errors 5)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -fmax-errors=${max_errors}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -ferror-limit=${max_errors}")
else()
	message(FATAL_ERROR "Unsupported compiler")
endif()

find_package(BISON REQUIRED)
find_package(RAGEL REQUIRED)

BISON_TARGET(ConfigParser config.yy ${CMAKE_CURRENT_BINARY_DIR}/parser.cpp COMPILE_FLAGS "-t -L c++")
RAGEL_TARGET(ConfigTokens config.rl ${CMAKE_CURRENT_BINARY_DIR}/lexer.cpp)
ADD_RAGEL_BISON_DEPENDENCY(ConfigTokens ConfigParser)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCES Compiler.cpp Polygon.cpp Simplex.cpp exprwrap.cpp main.cpp)

add_executable(main ${BISON_ConfigParser_OUTPUTS} ${RAGEL_ConfigTokens_OUTPUTS} ${SOURCES})

