
##
# CMake script for large project, with one executable
# and a library of objects
#
# The structure of the directory is assumed to be:
# ./source: containing implementations and one file with "int main()"
# ./include: containing all class declarations
# ./gtests: containing gtest unit tests
##

# Set the name of the project and target
# If your application follows the structure above, you don't need to
# specify anything else.
SET(_main_target praktikum)
SET(_main "")
#SET(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall")
add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fexceptions" "-Wfloat-equal" "-std=c++17"
  "$<$<CONFIG:DEBUG>:-O0;-g3;-ggdb>"
)
############################################################
# Normally you shouldn't need to change anything below.
############################################################
# Declare all source files the target consists of:
file(GLOB _files source/*cc include/*.h)

# Don't compile the main file into the library
IF(NOT "${_main}" STREQUAL "")
LIST(REMOVE_ITEM _files
${CMAKE_CURRENT_SOURCE_DIR}/${_main}
)
ENDIF()

INCLUDE_DIRECTORIES(include)

CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(${_main_target})

find_package (Eigen3 3.4 REQUIRED)
INCLUDE_DIRECTORIES( EIGEN3_INCLUDE_DIR )

#add_executable (example example.cpp)
#target_link_libraries (example Eigen3::Eigen)


OPTION(ENABLE_GOOGLE_TESTING
"Enable googletest style tests in this application." ON)

# We add one library and one target for  both Release and Debug
# mode, this will generate both Release and Debug programs for you.
# The debug library and program are postfixed with ".g"
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")

FOREACH(_build_type ${_d2_build_types})
  # Postfix to use everywhere
  SET(_p "${${_build_type}_postfix}")

    STRING(TOUPPER "${_build_type}" _BUILD_TYPE)

    # Setup library
    SET(_lib "${_main_target}${_p}")
    MESSAGE("-- Configuring library ${_lib}")

    ADD_LIBRARY(${_lib} SHARED ${_files})
    TARGET_LINK_LIBRARIES(${_lib}  Eigen3::Eigen )
    SET(TEST_LIBRARIES_${_BUILD_TYPE} ${_lib})

ENDFOREACH()

# Google test style testing
IF(ENABLE_GOOGLE_TESTING)
  find_package(GTest)
  if (GTest_FOUND)
      ADD_SUBDIRECTORY(gtests)
  else ()
      MESSAGE(WARNING "GTest library not found. Tests will not be compiled!")
  endif ()
ENDIF()

# Doxygen
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
    MESSAGE("-- Configured Doxygen")
  endif(DOXYGEN_FOUND)
endif()

add_custom_target(indent ./scripts/indent WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENABLE_TESTING()
