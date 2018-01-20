#
# Makefile for Project present like RTree
# Created by Vladimir Nedved 2017
# Apache License 2.0
#

BUILD_TYPE=release

CC=g++
NV=nvcc

PREPROC=-D CALC_POINT -D USE_CUDA

GDBFLAG=-ggdb
CLDFLAGS = -lm -lpthread -lcudart -std=c++0x 
TMP_CFLAGS= -c -Wall -march=native -pipe $(PREPROC) $(CLDFLAGS)

NVGDBFLAGS=--debug --device-debug
NVFLAGS=--gpu-architecture=compute_50 --gpu-code=compute_50,sm_50 $(PREPROC) --library=cudart
NVCLDFLAGS= -c -lm -lpthread -lcudart --std=c++11 $(PREPROC)

ifeq ($(BUILD_TYPE), )
	BUILD_TYPE=debug
endif

ifeq ($(BUILD_TYPE), release)
	GDBFLAG=
	NVGDBFLAGS=
	CFLAGS= -Ofast $(TMP_CFLAGS)
	NFLAGS= -O3 $(NVCLDFLAGS)
else
	CFLAGS= -Og $(TMP_CFLAGS)
	NFLAGS= -O0 $(NVCLDFLAGS)
endif

CXXFLAGS="$(CFLAGS)"

all: search2_cuda.o first.o search.o search2.o test.o
	@echo "=========== $(BUILD_TYPE) =============="
	$(CC) search2_cuda.o first.o search.o search2.o test.o -o ala_rtree $(CLDFLAGS) -L/usr/local/cuda-9.1/lib64
	
search2_cuda.o: search2.cu
	@echo "START: $(NV) $(NVGDBFLAGS) $(NVFLAGS) search2.cu"
	$(NV) $(NVGDBFLAGS) $(NFLAGS) $(NVFLAGS) search2.cu -o search2_cuda.o

first.o: first.cpp
	$(CC) $(GDBFLAG) $(CFLAGS) first.cpp
	#$(NV) $(NVGDBFLAGS) $(NFLAGS) first.cpp

search.o: search.cpp
	$(CC) $(GDBFLAG) $(CFLAGS) search.cpp
	
search2.o: search2.cpp
	$(CC) $(GDBFLAG) $(CFLAGS) search2.cpp
	
test.o: test.cpp
	$(CC) $(GDBFLAG) $(CFLAGS) test.cpp

clean:
	rm -rf *.o ala_rtree