CC=g++
CFLAGS=-c -Wall 
LDFLAGS=
SOURCES=main.cpp bch.cpp sourceCode/galois/GaloisFieldPolynomial.cpp sourceCode/galois/GaloisFieldElement.cpp sourceCode/galois/GaloisField.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@