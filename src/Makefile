CC=g++
LDFLAGS=-std=c++11 -O3 -lm
SOURCES=partitioner.cpp main.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=FM_Partitioner
INCLUDES=cell.h net.h partitioner.h

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o:  %.c  ${INCLUDES}
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)
