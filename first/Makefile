#
# Makefile for Project present like RTree
# Created by Vladimir Nedved 2017
# Apache License 2.0
#

CC=g++

GDBFLAG=-ggdb
CLFLAGS = -lm -lpthread
CFLAGS= -c -Wall

all: first.o search.o
	$(CC) first.o search.o -o ala_rtree $(CLFLAGS)

first.o: first.cpp
	$(CC) $(GDBFLAG) $(CFLAGS) $(CLFLAGS) first.cpp

search.o: search.cpp
	$(CC) $(GDBFLAG) $(CFLAGS) $(CLFLAGS) search.cpp

clean:
	rm -rf *.o ala_rtree