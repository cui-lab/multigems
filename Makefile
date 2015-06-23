CC=g++
CFLAGS=-c -O3 -Wall -Wno-sign-compare -std=c++0x -fopenmp
LDFLAGS= -fopenmp
SOURCES=gems.cpp core_functions.cpp multi_seq_obj.cpp seq_obj.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=multigems

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
.PHONY:clean
clean:
	rm -rf multigems
	rm -rf *.o
# DO NOT DELETE
# Debug mode
#
#CFLAGS=-c -Wall -Wno-sign-compare -pg -lc_p -std=c++0x -fopenmp
#LDFLAGS= -pg -fopenmp

