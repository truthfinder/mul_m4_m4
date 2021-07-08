@echo off

taskkill /IM PerfWatson2.exe /FI "STATUS eq RUNNING" /F

rd /S /Q build64
mkdir build64 && cd build64

rem mkdir build64 & pushd build64

rem "Visual Studio 15 2017 Wun64"
rem -A x64|Win32 -B "build32|build64" -S \path_to_source\
rem --build build64 --config Release
cmake .. -G "Visual Studio 16 2019" -A x64
cmake --build . --config Release --clean-first

rem popd
