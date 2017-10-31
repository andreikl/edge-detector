# CMake entry point
cmake_minimum_required (VERSION 2.6)
project (edge-detector)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

message( "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}" )
message( "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}" )
message( "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}" )
message( "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}" )

if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory." )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
    message( "Your Source Directory contains spaces." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
    message( "Your Build Directory contains spaces." )
endif()

# Compile external dependencies 
#add_subdirectory (external)

if(POLICY CMP0026)
    cmake_policy(SET CMP0026 OLD)
endif()

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/rpavlik-cmake-modules-fe2273")
#include(CreateLaunchers)
#include(MSVCMultipleProcessCompile) # /MP


include_directories(
     Streams/source
)

#set(ALL_LIBS
#    ${OPENGL_LIBRARY}
#    glfw
#    GLEW_1130
#)

#add_definitions(
#    -DTW_STATIC
#    -DTW_NO_LIB_PRAGMA
#    -DTW_NO_DIRECT3D
#    -DGLEW_STATIC
#    -D_CRT_SECURE_NO_WARNINGS
#    ${Qt5Widgets_DEFINITIONS}
#)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

add_executable(edge-detector
    src/utils.cpp
    src/main.cpp
)

target_link_libraries(edge-detector
    ${ALL_LIBS}
)

#set_target_properties(ogl PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
#create_target_launcher(ogl WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")

#SOURCE_GROUP(shaders REGULAR_EXPRESSION ".*/.*shader$" )

#add_custom_command(TARGET ogl POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different # which executes "cmake - E copy_if_different..."
#        "${_qt5Widgets_install_prefix}/bin/Qt5Widgetsd.dll"
#        $<TARGET_FILE_DIR:ogl>)

#add_custom_command(TARGET ogl POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different # which executes "cmake - E copy_if_different..."
#        "${_qt5Widgets_install_prefix}/bin/Qt5Cored.dll"
#        $<TARGET_FILE_DIR:ogl>)

#add_custom_command(TARGET ogl POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different # which executes "cmake - E copy_if_different..."
#        "${_qt5Widgets_install_prefix}/bin/Qt5Guid.dll"
#        $<TARGET_FILE_DIR:ogl>)
