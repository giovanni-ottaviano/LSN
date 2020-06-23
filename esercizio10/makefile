CXX=g++
MPICXX=mpicxx
CXXFLAGS=-std=c++11 -Wall -O

RANDOBJ=random.o
RANDDIR=Parallel_Generator
DEPS=$(RANDDIR)/random.h
GOAL=$(MAKECMDGOALS)
OBJECTS=$(GOAL).o $(RANDOBJ)


%.o: %.cpp $(DEPS)
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

$(RANDOBJ): $(RANDDIR)/random.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(GOAL): $(OBJECTS) $(DEPS)
	$(MPICXX) $(CXXFLAGS) $(OBJECTS) -o $(GOAL)
	rm $(OBJECTS)
