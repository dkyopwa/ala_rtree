#
# Makefile for Project present like RTree
# Created by Vladimir Nedved 2017
# Apache License 2.0
#

BUILD_TYPE=release

CC=g++

GDBFLAG=-ggdb
CLDFLAGS = -lm -lpthread
TMP_CFLAGS= -c -Wall -march=native -pipe $(CLDFLAGS)

ifeq ($(BUILD_TYPE), )
	BUILD_TYPE=debug
endif

ifeq ($(BUILD_TYPE), release)
	GDBFLAG=
	CFLAGS= -Ofast $(TMP_CFLAGS)
else
	CFLAGS= -Og $(TMP_CFLAGS)
endif

CXXFLAGS="$(CFLAGS)"

all: first.o search.o
	@echo "=========== $(BUILD_TYPE) =============="
	$(CC) first.o search.o -o ala_rtree $(CLDFLAGS)

first.o: first.cpp
	$(CC) $(GDBFLAG) $(CFLAGS) first.cpp

search.o: search.cpp
	$(CC) $(GDBFLAG) $(CFLAGS) search.cpp

clean:
	rm -rf *.o ala_rtree