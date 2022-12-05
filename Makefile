#!/bin/bash
OS_NAME := $(shell uname -s | tr A-Z a-z)


all: test 

os:
	@echo $(OS_NAME)

test: test.cpp
ifeq ($(OS_NAME),darwin)
	@echo $(OS_NAME)
	clang++ -std=c++17 -O3 test.cpp -o test.x -framework accelerate
else 
ifeq ($(OS_NAME),linux)
	@echo $(OS_NAME)
	g++ -std=c++17 -O3 test.cpp -o test.x
else
	@echo $(OS_NAME)
	g++ -std=c++17 -O3 test.cpp -o test.x
endif
endif

clean: 
	rm -f *.o *.x