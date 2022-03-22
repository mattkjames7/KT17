@echo off

echo Compiling model

g++ -lm -fopenmp -fPIC -std=c++17 -c dipole.cc -o dipole.o
g++ -lm -fopenmp -fPIC -std=c++17 -c disk.cc -o disk.o
g++ -lm -fopenmp -fPIC -std=c++17 -c qhsheet.cc -o qhsheet.o
g++ -lm -fopenmp -fPIC -std=c++17 -c shield.cc -o shield.o
g++ -lm -fopenmp -fPIC -std=c++17 -c kt14.cc -o kt14.o
g++ -lm -fopenmp -fPIC -std=c++17 -c kt17.cc -o kt17.o
g++ -lm -fopenmp -fPIC -std=c++17 -c ktmodel.cc -o ktmodel.o
g++ -lm -fopenmp -fPIC -std=c++17 -c model.cc -o model.o

exit /b 0

:CompileError
echo Compilation error
exit /b 8
