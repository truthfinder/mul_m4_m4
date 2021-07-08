cmake_minimum_required(VERSION 3.18)
project(mul_m4_m4)

file(GLOB SOURCES "*.h" "*.cpp")
source_group("" FILES ${SOURCES})

add_definitions(-DIACA_MARKS_OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /c /std:c++17 /GL /arch:AVX /EHsc /D_WIN64")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")

add_executable(prog-mtx WIN32 ${SOURCES})

#set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} /SUBSYSTEM:CONSOLE /LTCG")
set_target_properties(prog-mtx PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE /LTCG")

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT prog-mtx)
#CMAKE_VS_PLATFORM_NAME_DEFAULT
