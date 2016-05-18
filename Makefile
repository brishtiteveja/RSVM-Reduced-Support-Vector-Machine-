CC = gcc
CXXC = g++
CFLAGS = -Wall -O3 -g

all: svm-train svm-predict svm-scale

#solvebqp1.o: solvebqp1.c svm.h
#	$(CC) $(CFLAGS) -c solvebqp1.c

svm-predict: svm-predict.c svm.o solvebqp.o
	$(CXXC) $(CFLAGS) svm-predict.c svm.o solvebqp.o dtron/tron.a f2c/f2c.a -o svm-predict -lm
svm-train: svm-train.c svm.o solvebqp.o
	$(CXXC) $(CFLAGS) svm-train.c svm.o solvebqp.o dtron/tron.a f2c/f2c.a -o svm-train -lm
svm-scale: svm-scale.c
	$(CXXC) $(CFLAGS) svm-scale.c -o svm-scale
#	$(CC) $(CFLAGS) svm-train.c svm.o solvebqp.o dtron/tron.a f2c/ibm.a -o svm-train -lm
svm.o: svm.cpp svm.h
	$(CXXC) $(CFLAGS) -c svm.cpp
solvebqp.o: solvebqp.c svm.h
	$(CC) $(CFLAGS) -c solvebqp.c
commit:
	LOGNAME=adm cvs commit
clean:
	rm -f *~ *.o svm-train svm-predict svm-scale
