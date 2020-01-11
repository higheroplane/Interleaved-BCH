CC=g++
CFLAGS=-c -std=c++11
LDFLAGS=
SOURCES=main.cpp bch.cpp linalg.cpp stoch.cpp  \
sourceCode/galois/GaloisFieldPolynomial.cpp sourceCode/galois/GaloisFieldElement.cpp sourceCode/galois/GaloisField.cpp \
Bigint/karatsuba.cpp BigInt/bigint.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o