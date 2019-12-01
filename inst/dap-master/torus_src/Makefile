GSLV := $(shell pkg-config --libs gsl)
GSLI := $(shell pkg-config --cflags gsl)
BOOSTL := -L$(shell echo "${LIBRARY_PATH}" | tr ':' '\n' | grep 'boost')
BOOSTI := -I$(shell echo "${CPATH}" | tr ':' '\n' | grep 'boost')

main: main.o controller.o logistic.o
	$(CXX)  $(CXXFLAGS)  -g main.o controller.o logistic.o  -lm $(GSLV) $(BOOSTL) -lboost_iostreams -lz  -o torus
static: main.o controller.o logistic.o
	$(CXX)  $(CXXFLAGS) $(GSLI) $(BOOSTI) -g main.o controller.o logistic.o  -lm $(GSLV) $(BOOSTL) -lboost_iostreams -lz -static -o torus.static
main.o: main.cc 
	$(CXX) $(CXXFLAGS) -g $(GSLI) $(BOOSTI) -c  main.cc
controller.o: controller.cc classdef.h
	$(CXX) $(CXXFLAGS)  -g $(GSLI) $(BOOSTI) -c  controller.cc
logistic.o: logistic.cc logistic.h
	$(CXX) $(CXXFLAGS) -g $(GSLI) $(BOOSTI) -c logistic.cc -fpermissive
clean:
	rm *.o torus	
