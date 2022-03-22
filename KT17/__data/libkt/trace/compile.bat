@echo off

echo Compiling trace

g++ -lm -fopenmp -fPIC -std=c++17 -c conttrace.cc -o conttrace.o
g++ -lm -fopenmp -fPIC -std=c++17 -c interptraceclosestpos.cc -o interptraceclosestpos.o
g++ -lm -fopenmp -fPIC -std=c++17 -c latlt.cc -o latlt.o
g++ -lm -fopenmp -fPIC -std=c++17 -c fieldlinernorm.cc -o fieldlinernorm.o
g++ -lm -fopenmp -fPIC -std=c++17 -c trace.cc -o trace.o
	
exit /b 0

:CompileError
echo Compilation error
exit /b 8
