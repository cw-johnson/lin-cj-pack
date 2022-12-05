#!/bin/bash
all: test

test: test.cpp
	clang++ -std=c++17 -O3 test.cpp -o test.x -framework accelerate

clean: 
	rm -f *.o *.x