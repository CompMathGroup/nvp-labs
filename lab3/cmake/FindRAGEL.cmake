#.rst:
# FindRAGEL
# --------
#
# Find ragel executable and provides a macro to generate custom build rules
#
#
#
# The module defines the following variables:
#
# ::
#
#   RAGEL_FOUND - true is ragel executable is found
#   RAGEL_EXECUTABLE - the path to the ragel executable
#   RAGEL_VERSION - the version of ragel
#
#
#
# The minimum required version of ragel can be specified using the
# standard syntax, e.g.  find_package(RAGEL 6.8)
#
#
#
# If ragel is found on the system, the module provides the macro:
#
# ::
#
#   RAGEL_TARGET(Name RagelInput RagelOutput [COMPILE_FLAGS <string>])
#
# which creates a custom command to generate the <RagelOutput> file from
# the <RagelInput> file.  If COMPILE_FLAGS option is specified, the next
# parameter is added to the ragel command line.  Name is an alias used to
# get details of this custom command.  Indeed the macro defines the
# following variables:
#
# ::
#
#   RAGEL_${Name}_DEFINED - true is the macro ran successfully
#   RAGEL_${Name}_OUTPUTS - the source file generated by the custom rule, an
#   alias for RagelOutput
#   RAGEL_${Name}_INPUT - the ragel source file, an alias for ${RagelInput}
#
#
#
# Ragel scanners may use tokens defined by Bison: the code generated
# by Ragel depends of the header generated by Bison.  This module also
# defines a macro:
#
# ::
#
#   ADD_RAGEL_BISON_DEPENDENCY(RagelTarget BisonTarget)
#
# which adds the required dependency between a scanner and a parser
# where <RagelTarget> and <BisonTarget> are the first parameters of
# respectively RAGEL_TARGET and BISON_TARGET macros.
#
# ::
#
#   ====================================================================
#   Example:
#
#
#
# ::
#
#    find_package(BISON)
#    find_package(RAGEL)
#
#
#
# ::
#
#    BISON_TARGET(MyParser parser.y ${CMAKE_CURRENT_BINARY_DIR}/parser.cpp)
#    RAGEL_TARGET(MyScanner lexer.rl  ${CMAKE_CURRENT_BINARY_DIR}/lexer.cpp)
#    ADD_RAGEL_BISON_DEPENDENCY(MyScanner MyParser)
#
#
#
# ::
#
#    include_directories(${CMAKE_CURRENT_BINARY_DIR})
#    add_executable(Foo
#       Foo.cc
#       ${BISON_MyParser_OUTPUTS}
#       ${RAGEL_MyScanner_OUTPUTS}
#    )
#   ====================================================================

#=============================================================================
# Copyright 2009 Kitware, Inc.
# Copyright 2006 Tristan Carel
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

find_program(RAGEL_EXECUTABLE NAMES ragel DOC "path to the ragel executable")
mark_as_advanced(RAGEL_EXECUTABLE)

if(RAGEL_EXECUTABLE)

  execute_process(COMMAND ${RAGEL_EXECUTABLE} --version
    OUTPUT_VARIABLE RAGEL_version_output
    ERROR_VARIABLE RAGEL_version_error
    RESULT_VARIABLE RAGEL_version_result
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT ${RAGEL_version_result} EQUAL 0)
    if(RAGEL_FIND_REQUIRED)
      message(SEND_ERROR "Command \"${RAGEL_EXECUTABLE} --version\" failed with output:\n${RAGEL_version_output}\n${RAGEL_version_error}")
    else()
      message("Command \"${RAGEL_EXECUTABLE} --version\" failed with output:\n${RAGEL_version_output}\n${RAGEL_version_error}\nRAGEL_VERSION will not be available")
    endif()
  else()
    string(REGEX REPLACE "^.*version ([0-9]+[^ ]*).*$" "\\1"
      RAGEL_VERSION "${RAGEL_version_output}")
    unset(RAGEL_EXE_EXT)
    unset(RAGEL_EXE_NAME_WE)
  endif()

  #============================================================
  # RAGEL_TARGET (public macro)
  #============================================================
  #
  macro(RAGEL_TARGET Name Input Output)
    set(RAGEL_TARGET_usage "RAGEL_TARGET(<Name> <Input> <Output> [COMPILE_FLAGS <string>]")
    if(${ARGC} GREATER 3)
      if(${ARGC} EQUAL 5)
        if("${ARGV3}" STREQUAL "COMPILE_FLAGS")
          set(RAGEL_EXECUTABLE_opts  "${ARGV4}")
          separate_arguments(RAGEL_EXECUTABLE_opts)
        else()
          message(SEND_ERROR ${RAGEL_TARGET_usage})
        endif()
      else()
        message(SEND_ERROR ${RAGEL_TARGET_usage})
      endif()
    endif()

    add_custom_command(OUTPUT ${Output}
      COMMAND ${RAGEL_EXECUTABLE}
      ARGS ${RAGEL_EXECUTABLE_opts} -o${Output} ${Input}
      DEPENDS ${Input}
      COMMENT "[RAGEL][${Name}] Building scanner with ragel ${RAGEL_VERSION}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    set(RAGEL_${Name}_DEFINED TRUE)
    set(RAGEL_${Name}_OUTPUTS ${Output})
    set(RAGEL_${Name}_INPUT ${Input})
    set(RAGEL_${Name}_COMPILE_FLAGS ${RAGEL_EXECUTABLE_opts})
  endmacro()
  #============================================================


  #============================================================
  # ADD_RAGEL_BISON_DEPENDENCY (public macro)
  #============================================================
  #
  macro(ADD_RAGEL_BISON_DEPENDENCY RagelTarget BisonTarget)

    if(NOT RAGEL_${RagelTarget}_OUTPUTS)
      message(SEND_ERROR "Ragel target `${RagelTarget}' does not exists.")
    endif()

    if(NOT BISON_${BisonTarget}_OUTPUT_HEADER)
      message(SEND_ERROR "Bison target `${BisonTarget}' does not exists.")
    endif()

    set_source_files_properties(${RAGEL_${RagelTarget}_OUTPUTS}
      PROPERTIES OBJECT_DEPENDS ${BISON_${BisonTarget}_OUTPUT_HEADER})
  endmacro()
  #============================================================

endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(RAGEL REQUIRED_VARS RAGEL_EXECUTABLE
                                       VERSION_VAR RAGEL_VERSION)

# FindRAGEL.cmake ends here
