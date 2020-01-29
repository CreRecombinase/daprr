GSL_LIBS := $(shell gsl-config --libs)
GSL_INC := $(shell gsl-config --cflags)
CXX := $(shell R CMD config CXX)
CXX11STD := $(shell R CMD config CXX11STD)
main: main.o controller.o logistic.o 
	$(CXX)  $(CXX11STD) -O3 main.o controller.o logistic.o  $(GSL_LIBS)   -lboost_iostreams -lz  -o torus
static: main.o controller.o logistic.o
	$(CXX)  $(CXX11STD) -O3 main.o controller.o logistic.o   $(GSL_LIBS) -lboost_iostreams -lz -static -o torus.static
main.o: main.cc 
	$(CXX) $(CXX11STD) -c  main.cc $(GSL_INC)
controller.o: controller.cc classdef.h
	$(CXX) $(CXX11STD) -fopenmp  -c  controller.cc $(GSL_INC)
logistic.o: logistic.cc logistic.h 
	$(CXX) $(CXX11STD) -c logistic.cc -fpermissive $(GSL_INC)
clean:
	rm *.o torus	
