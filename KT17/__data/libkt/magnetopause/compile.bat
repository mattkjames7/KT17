@echo off

echo Compiling magnetopause

g++ -lm -fopenmp -fPIC -std=c++17 -c magnetopause.cc -o magnetopause.o
g++ -lm -fopenmp -fPIC -std=c++17 -c withinmp.cc -o withinmp.o

exit /b 0

:CompileError
echo Compilation error
exit /b 8
