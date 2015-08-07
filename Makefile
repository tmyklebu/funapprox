LDFLAGS=-lgmp -lgmpxx
LIBRARIES=QSex090408/lib/QSopt_ex.a
CXXFLAGS=-O3 -Wall -Wextra -g -std=gnu++0x -march=native -mtune=native -IQSex090408/include -Wno-unused -Wno-unused-parameter

.PHONY: all
all: fa

%: %.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LIBRARIES)
