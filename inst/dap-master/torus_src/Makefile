main: main.o controller.o logistic.o
	$(CXX)  $(CXXFLAGS) main.o controller.o logistic.o  -lm -lgsl -lgslcblas -lboost_iostreams -lz  -o torus
static: main.o controller.o logistic.o
	$(CXX)  $(CXXFLAGS) main.o controller.o logistic.o  -lm -lgsl -lgslcblas -lboost_iostreams -lz -static -o torus.static
main.o: main.cc 
	$(CXX) $(CXXFLAGS) -c  main.cc
controller.o: controller.cc classdef.h
	$(CXX) $(CXXFLAGS)  -c  controller.cc
logistic.o: logistic.cc logistic.h
	$(CXX) $(CXXFLAGS) -c logistic.cc -fpermissive
clean:
	rm *.o torus	
