CC = g++
CFLAGS = -O3 -march=native -fopenmp -funroll-loops -mfpmath=sse -DNDEBUG -DEIGEN_NO_DEBUG

sample: main.o FermionicFockState.o
	$(CC) $(CFLAGS) -o main.exe main.o FermionicFockState.o
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp -I "C:\Users\matte\C++_libraries\eigen-3.4.0\eigen-3.4.0"
FermionicFockState.o: FermionicFockState.cpp FermionicFockState.h
	$(CC) $(CFLAGS) -c FermionicFockState.cpp -I "C:\Users\matte\Documents\Study_C++\SampleProject" 
