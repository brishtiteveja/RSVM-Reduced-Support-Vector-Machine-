CC = gcc
CXXC = g++
CFLAGS = -Wall -g -O3

all: bsvm-train bsvm-predict svm-scale

bsvm-predict: svm-predict.c bsvm.o solvebqp.o tron
	$(CXXC) $(CFLAGS) svm-predict.c bsvm.o solvebqp.o \
dtron/tron.a f2c/f2c.a -o bsvm-predict -lm
bsvm-train: svm-train.c bsvm.o solvebqp.o tron 
	$(CXXC) $(CFLAGS) svm-train.c bsvm.o solvebqp.o \
dtron/tron.a f2c/f2c.a -o bsvm-train -lm
svm-scale: svm-scale.c
	$(CXXC) $(CFLAGS) svm-scale.c -o svm-scale
bsvm.o: bsvm.cpp svm.h
	$(CXXC) $(CFLAGS) -c bsvm.cpp
solvebqp.o: solvebqp.c svm.h
	$(CC) $(CFLAGS) -c solvebqp.c
tron: svm.h
	cd dtron; make OPTFLAGS='$(CFLAGS)' CC='$(CC)';
	cd f2c; make OPTFLAGS='$(CFLAGS)' CC='$(CC)';
clean:
	cd f2c;		make clean
	cd dtron;	make clean
	rm -f *~ bsvm.o solvebqp.o bsvm-train bsvm-predict svm-scale
