CC=g++ -march=native -O3
CFLAGS=-c -I. -std=c++11 -Wfatal-errors

all: dsample

dsample: dsample.o ddecom.o 
	$(CC) dsample.o ddecom.o -o dsample
	rm *.o

dsample.o: dsample.cc
	$(CC) $(CFLAGS) dsample.cc -o dsample.o

ddecom.o: ddecom.cc
	$(CC) $(CFLAGS) ddecom.cc -o ddecom.o
