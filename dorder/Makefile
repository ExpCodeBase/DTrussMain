CC=g++ -march=native -O3
CFLAGS=-c -I. -std=c++0x -Wfatal-errors

all: dm

dm: dtest.o dorder.o
	$(CC) dtest.o dorder.o -o dm
	rm *.o

dtest.o: dtest.cc
	$(CC) $(CFLAGS) dtest.cc -o dtest.o

dorder.o: dorder.cc
	$(CC) $(CFLAGS) dorder.cc -o dorder.o
