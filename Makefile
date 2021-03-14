#!bin/make

compile: src/main.cpp src/variational_bayes.cpp src/SigTracer.h src/Sampler.h
	g++ -O2 src/main.cpp -o bin/ST -std=c++11

