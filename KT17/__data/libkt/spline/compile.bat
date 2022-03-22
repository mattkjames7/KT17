@echo off

echo Compiling spline

g++ -lm -fopenmp -fPIC -std=c++17 -c spline.cc -o spline.o
g++ -lm -fopenmp -fPIC -std=c++17 -c libspline.cc -o libspline.o

exit /b 0

:CompileError
echo Compilation error
exit /b 8
