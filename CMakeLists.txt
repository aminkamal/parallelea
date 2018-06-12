cmake_minimum_required (VERSION 3.0 FATAL_ERROR)

# Change this section to correspond with your project info.
project(smart-o2cpp LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE OFF)

get_filename_component (PROJECT_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src" REALPATH)
get_filename_component (HEADERS "include" REALPATH)
configure_file ("${HEADERS}/smarto2cpp.h.in" "${HEADERS}/smarto2cpp.h")
include_directories (${HEADERS})

set (CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem")

set (PROJECT_SUMMARY                            "Strathclyde Mechanical and Aerospace Researt Toolbox for Optimisation and Optimal Control (SMART-O2CPP)")
set (${CMAKE_PROJECT_NAME}_VENDOR_NAME          "spaceART")
set (${CMAKE_PROJECT_NAME}_VENDOR_CONTACT       "spaceART@strathclyde.ac.uk")

set (PROJECT_PATH                               "${CMAKE_CURRENT_SOURCE_DIR}")
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}     "${PROJECT_PATH}/cmake/Modules")
set (INCLUDE_PATH                               "${PROJECT_PATH}/include")
if (NOT DOCS_PATH)
  set (DOCS_PATH                                "${PROJECT_PATH}/docs")
endif()
set (LIB_PATH                                   "${PROJECT_BINARY_DIR}/lib")
set (LIB_NAME                                   "smart-o2cpp")
set (TEST_PATH                                  "${PROJECT_BINARY_DIR}/test")


option(BUILD_PARALLEL_EXAMPLE "Build pool and island-based evolutionary algorithms example" ON)
option(ENABLE_MULTITHREADED_POOL_IM "Parallelise pool and island-based evolutionary algorithms" ON)

SET (MANDATORY_LIBRARIES "")


if (ENABLE_MULTITHREADED_POOL_IM)
    add_definitions(-DUSEMULTITHREADS)
endif ()

if (BUILD_PARALLEL_EXAMPLE)
    add_executable(parallel examples/parallel_example.cpp)
    target_link_libraries(parallel ${LIB_NAME} ${MANDATORY_LIBRARIES})
endif ()

if (ENABLE_MULTITHREADED_POOL_IM)
    find_package(Threads REQUIRED)
endif ()

get_filename_component (HEADERS "include" REALPATH)
include_directories (${HEADERS})

file (GLOB_RECURSE SOURCES
        src/Parallel/*.cpp)
file (GLOB_RECURSE TEST_SOURCES "src/test/*.cpp")

#file (GLOB_RECURSE LIB_HEADERS "include/*.h")
file (GLOB_RECURSE LIB_HEADERS "include/Parallel/*.h")
add_custom_target (headers SOURCES ${LIB_HEADERS})


set (CMAKE_BUILD_TYPE Debug)

# Set CMake build-type. If it not supplied by the user, the default built type is "Release".
#if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
#    message(STATUS "WARNING: building release version!")
#    set (CMAKE_BUILD_TYPE Release)
#elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
#    message (STATUS "WARNING: building debug version!")
#endif()

# Set platform-specific compiler flags.
if (WIN32)
  if (MSVC)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    # We need to force MSVC to use static runtime:
    # (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace).
    foreach (flag_var
      CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      # Find all dynamic runtime (MD) references and replace with static (MT)
      if (${flag_var} MATCHES "/MD")
        string (REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif()
    endforeach()
  else (MSVC)
    set (CMAKE_LEGACY_CYGWIN_WIN32 0)
    set (CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Wnon-virtual-dtor")
  endif()
else()
  set (CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Wnon-virtual-dtor")
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
    #set (CMAKE_CXX_FLAGS_DEBUG   "-fprofile-arcs -ftest-coverage -O0 -g3")
    set (CMAKE_CXX_FLAGS_DEBUG   "-g3 -O3")
    set (CMAKE_CXX_FLAGS_RELEASE "-O3")
    set (CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}")
    if(CMAKE_USE_PTHREADS_INIT)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    endif ()
endif()

configure_file("${INCLUDE_PATH}/config.h.in" "${INCLUDE_PATH}/config.h")

if (BUILD_STATIC) 
   add_library (${LIB_NAME} STATIC ${SOURCES})
else()
   add_library (${LIB_NAME} SHARED ${SOURCES})
endif()

set_target_properties (${LIB_NAME} PROPERTIES
	  VERSION ${PROJECT_VERSION}
	  SOVERSION ${PROJECT_VERSION}
)

#add_executable (main src/main.cpp)
#target_link_libraries (main ${MANDATORY_LIBRARIES})
#add_dependencies(main ${LIB_NAME})

if (BUILD_DOCS)
  find_package (Doxygen)

  if (NOT DOXYGEN_FOUND)
    message (FATAL_ERROR
            "Doxygen is needed to build the documentation. Please make sure it is installed.")
  endif()

  configure_file (${PROJECT_PATH}/Doxyfile.in
                 ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  # Add a custom target to run Doxygen whenever the project is built.
  # If you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the command below.
  add_custom_target (docs ALL COMMAND ${DOXYGEN_EXECUTABLE}
                    ${PROJECT_BINARY_DIR}/Doxyfile
                    SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

#add_subdirectory ("examples")

# Install header files and library.
# Destination is set by CMAKE_INSTALL_PREFIX and defaults to usual locations, unless overridden by
# user.
install (DIRECTORY "${INCLUDE_PATH}/" DESTINATION ${CMAKE_PROJECT_NAME}/include FILES_MATCHING PATTERN "*.h")
install (TARGETS ${LIB_NAME} DESTINATION ${CMAKE_PROJECT_NAME}/lib)
install (TARGETS ${BIN_NAME} DESTINATION ${CMAKE_PROJECT_NAME}/bin)

# Set up packager.
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CMAKE_PROJECT_NAME}")
set (CPACK_PACKAGE_VENDOR              "${CMAKE_PROJECT_NAME}_VENDOR_NAME")
set (CPACK_PACKAGE_DESCRIPTION_FILE    "${PROJECT_PATH}/README.md")
set (CPACK_RESOURCE_FILE_LICENSE       "${PROJECT_PATH}/LICENSE.md")
set (CPACK_PACKAGE_VERSION_MAJOR       "${PROJECT_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR       "${PROJECT_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH       "${PROJECT_VERSION_PATCH}")
set (CPACK_PACKAGE_INSTALL_DIRECTORY   "${CMAKE_PROJECT_NAME}")
if (WIN32)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  set (CPACK_NSIS_DISPLAY_NAME         "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
  set (CPACK_NSIS_CONTACT              "${CMAKE_PROJECT_NAME}_VENDOR_CONTACT")
  set (CPACK_NSIS_MODIFY_PATH ON)
endif()

include (CPack)
