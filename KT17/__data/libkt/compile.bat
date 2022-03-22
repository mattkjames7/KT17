@echo off

where /q gcc
if %ERRORLEVEL% neq 0 (
	echo GCC wasn't found
	exit /b 6
)

WHERE /Q gfortran 
if %ERRORLEVEL% neq 0 (
	echo GFortran wasn't found 
	exit /b 7
)

echo Compiling libkt...

cd magnetopause
call compile.bat
cd ..

cd model
call compile.bat
cd ..

cd spline
call compile.bat
cd ..

cd trace
call compile.bat
cd ..


g++ -fPIC -c -lm -fopenmp -std=c++17 libkt.cc -o libkt.o
if %ERRORLEVEL% neq 0 (goto CompileError)

g++ -lm -fopenmp -fPIC -std=c++17 -shared -o libkt.dll *.o magnetopause/*.o model/*.o spline/*.o	trace/*.o
if %ERRORLEVEL% neq 0 (goto CompileError)


echo Done

del *.o
del magnetopause\*.o
del model\*.o
del spline\*.o
del trace\*.o
exit /b 0

:CompileError
echo Compilation error
exit /b 8
