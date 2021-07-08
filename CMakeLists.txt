cmake_minimum_required(VERSION 3.18)
project(prog-mtx)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#set_target_properties(prog-mtx PROPERTIES
#			CXX_STANDARD 17
#			CXX_STANDARD_REQUIRED ON
#			CXX_EXTENSIONS OFF)

file(GLOB SOURCES "*.h" "*.cpp")
source_group("" FILES ${SOURCES})

add_definitions(-DIACA_MARKS_OFF)

if(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /GL /arch:AVX /EHsc /D_WIN64")
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest")

	add_executable(${PROJECT_NAME} WIN32 ${SOURCES}) #WIN32->EXCLUDE_FROM_ALL

	set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE /LTCG")
else(WIN32)
	add_executable(${PROJECT_NAME} ${SOURCES})
endif(WIN32)

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT prog-mtx)
#CMAKE_VS_PLATFORM_NAME_DEFAULT
