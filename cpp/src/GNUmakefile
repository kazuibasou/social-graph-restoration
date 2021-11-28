CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

all:mkdir main clean

mkdir:
	mkdir -p ../../gen_graph
	mkdir -p ../bin

main: basic_function.o graph.o sampling.o estimation.o rewiring.o generation.o main.o
	$(CXX) $(CXXFLAGS) -o ../bin/main basic_function.o graph.o sampling.o estimation.o rewiring.o generation.o main.o -O3
main.o: main.cpp basic_function.h graph.h sampling.h estimation.h rewiring.h generation.h
	$(CXX) $(CXXFLAGS) -c main.cpp -O3
basic_function.o:basic_function.cpp basic_function.h
	$(CXX) $(CXXFLAGS) -c basic_function.cpp -O3
graph.o:graph.cpp graph.h basic_function.h
	$(CXX) $(CXXFLAGS) -c graph.cpp -O3
sampling.o:sampling.cpp basic_function.h graph.h sampling.h
	$(CXX) $(CXXFLAGS) -c sampling.cpp -O3
estimation.o:estimation.cpp basic_function.h graph.h estimation.h
	$(CXX) $(CXXFLAGS) -c estimation.cpp -O3
rewiring.o:rewiring.cpp basic_function.h graph.h rewiring.h
	$(CXX) $(CXXFLAGS) -c rewiring.cpp -O3
generation.o:generation.cpp basic_function.h graph.h estimation.h rewiring.h generation.h
	$(CXX) $(CXXFLAGS) -c generation.cpp -O3

clean:
	rm *.o