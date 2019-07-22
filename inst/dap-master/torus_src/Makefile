main: main.o controller.o logistic.o 
	g++  -std=c++11 -O3 main.o controller.o logistic.o  -lm -lgsl -lgslcblas -lboost_iostreams -lz  -o torus
static: main.o controller.o logistic.o
	g++  -std=c++11 -O3 main.o controller.o logistic.o  -lm -lgsl -lgslcblas -lboost_iostreams -lz -static -o torus.static
main.o: main.cc 
	g++ -std=c++11 -c  main.cc
controller.o: controller.cc classdef.h
	g++ -std=c++11  -c  controller.cc
logistic.o: logistic.cc logistic.h
	g++ -std=c++11 -c logistic.cc -fpermissive
clean:
	rm *.o torus	
