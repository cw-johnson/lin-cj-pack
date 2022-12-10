#!/bin/bash
OS_NAME := $(shell uname -s | tr A-Z a-z)


all: hw2 hw3 

os:
	@echo $(OS_NAME)


setup:
ifeq ($(OS_NAME),darwin)
	@echo $(OS_NAME)
	xcode-select --install
else 
ifeq ($(OS_NAME),linux)
	@echo $(OS_NAME)
	sudo apt-get install libblas-dev liblapack-dev
else
	@echo $(OS_NAME)
endif
endif

hw2: hw2.cpp
ifeq ($(OS_NAME),darwin)
	@echo $(OS_NAME)
	clang++ -std=c++20 -O3 hw2.cpp -o hw2.x -framework accelerate
else 
ifeq ($(OS_NAME),linux)
	@echo $(OS_NAME)
	g++ -std=c++20 -O3 hw2.cpp -o hw2.x -lblas
else
	@echo $(OS_NAME)
	g++ -std=c++20 -O3 hw2.cpp -o hw2.x -lblas
endif
endif

hw3: hw3.cpp
ifeq ($(OS_NAME),darwin)
	@echo $(OS_NAME)
	clang++ -std=c++20 -O3 hw3.cpp -o hw3.x -framework accelerate
else 
ifeq ($(OS_NAME),linux)
	@echo $(OS_NAME)
	g++ -std=c++20 -O3 hw3.cpp -o hw3.x -lblas
else
	@echo $(OS_NAME)
	g++ -std=c++20 -O3 hw3.cpp -o hw3.x -lblas
endif
endif

clean: 
	rm -f *.o *.x