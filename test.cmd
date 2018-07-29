@rem Enterprise, Professional или Community
@rem \Program Files (x86)\Microsoft Visual Studio\Version\Offering\Common7\Tools
@rem \Program Files\Microsoft Visual Studio\Version\Offering\Common7\Tools
@rem VsDevCmd.bat
@rem C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat
@rem if "%1"=="x32" call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars32.bat"
@rem if "%1"=="x64" call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat"
@rem if "%1"=="" call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat" -arch=x86

@echo off

set arch="x32"
set iaca=0
set iaca2=0
set defines="/D__SSE__ /D__AVX__"
set fname="prog-mtx"

set iarch="HSW"
rem set iarch="SKX"

for %%x in (%*) do (
if %%x==x32 (set arch="x32")
if %%x==x86 (set arch="x32")
if %%x==x64 (set arch="x64")
if %%x==avx (set defines=%defines% "/D__AVX__")
if %%x==fma (set defines=%defines% "/D__AVX__ /D__FMA__")
if %%x==avx2 (set defines=%defines% "/D__AVX__ /D__FMA__" "/D__AVX2__")
if %%x==avx512 (set defines=%defines% "/D__AVX__ /D__FMA__" "/D__AVX2__" "/D__AVX512__")
if %%x==iaca (set iaca=1)
if %%x==iaca3 (set iaca=1)
if %%x==iaca2 (set iaca2=1)
if %%x==ivb (set iarch="IVB")
if %%x==hsw (set iarch="HSW")
if %%x==bdw (set iarch="BDW")
if %%x==skl (set iarch="SKL")
if %%x==skx (set iarch="SKX")
)

@rem if not exist filename.txt exit /b 1

if %arch%=="x32" call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat" -arch=x86
if %arch%=="x64" call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat" -arch=amd64

del *.exe
del %fname%*.txt
del *.asm

if %iaca%==1 if %arch%=="x32" goto :do_iaca_x32
if %iaca%==1 if %arch%=="x64" goto :do_iaca_x64
if %iaca2%==1 if %arch%=="x32" goto :do_iaca2_x32
if %iaca2%==1 if %arch%=="x64" goto :do_iaca2_x64
if %iaca%==0 if %arch%=="x32" goto :do_x32
if %iaca%==0 if %arch%=="x64" goto :do_x64
goto :error

rem cl GR- /arch:AVX2 /Fatest.asm /FAsu /link /NOLOGO /MACHINE:X64

:do_x32
echo "compile x32"
cl /std:c++17 /O2 /GL /arch:AVX /EHsc /DIACA_MARKS_OFF %defines% /Fa%fname%.asm /FAsu %fname%.cpp /link /NOLOGO > %fname%32.txt
if %errorlevel% neq 0 (echo "cl error" && goto :error)
%fname%.exe
goto :end

:do_x64
echo "compile x64"
cl /std:c++17 /O2 /GL /arch:AVX /EHsc /D_WIN64 /DIACA_MARKS_OFF %defines% /Fa%fname%.asm /FAsu /Fe%fname%64.exe %fname%.cpp /link /NOLOGO /MACHINE:X64 /OPT:REF /OPT:ICF > %fname%64.txt
if %errorlevel% neq 0 (echo "cl error" && goto :error)
%fname%64.exe
goto :end

:do_iaca_x32
echo "iaca3 x32"
cl /c /std:c++17 /O2 /GL /arch:AVX /EHsc /DIACA3 %defines% /Fa%fname%.asm /FAsu %fname%.cpp > %fname%32.txt
if %errorlevel% neq 0 (echo "cl error" && goto :error)
..\iaca\iaca.exe -arch %iarch% %fname%.obj > %fname%32_iaca.txt
if %errorlevel% neq 0 (echo "iaca error" && goto :error)
goto :end

:do_iaca_x64
echo "iaca3 x64"
cl /c /std:c++17 /O2 /GL /arch:AVX /EHsc /DIACA3 /D_WIN64 %defines% /Fa%fname%.asm /FAsu %fname%.cpp > %fname%64.txt
if %errorlevel% neq 0 (echo "cl error" && goto :error)
..\iaca\iaca.exe -arch %iarch% %fname%.obj > %fname%64_iaca.txt
if %errorlevel% neq 0 (echo "iaca error" && goto :error)
goto :end

:do_iaca2_x32
echo "iaca2 x32"
cl /c /std:c++17 /O2 /GL /arch:AVX /EHsc /DIACA2 %defines% /Fa%fname%.asm /FAsu %fname%.cpp > %fname%32.txt
if %errorlevel% neq 0 (echo "cl error" && goto :error)
..\iaca2\iaca.exe -32 -arch %iarch% %fname%.obj > %fname%32_iaca.txt
if %errorlevel% neq 0 (echo "iaca error" && goto :error)
goto :end

:do_iaca2_x64
echo "iaca3 x64"
cl /c /std:c++17 /O2 /GL /arch:AVX /EHsc /DIACA2 /D_WIN64 %defines% /Fa%fname%.asm /FAsu %fname%.cpp > %fname%64.txt
if %errorlevel% neq 0 (echo "cl error" && goto :error)
..\iaca2\iaca.exe -64 -arch %iarch% %fname%.obj > %fname%64_iaca.txt
if %errorlevel% neq 0 (echo "iaca error" && goto :error)
goto :end

:end
del *.obj
goto :eof

:error
echo "script failed"
exit /b 1

