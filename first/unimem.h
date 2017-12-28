/*
	Project present like RTree
	Created by Vladimir Nedved 2017
	Apache License 2.0
	file unimem.h
*/
#pragma once
#include "stdlib.h"
#include "malloc.h"

#ifdef _WIN
inline void* aligned_alloc(size_t alignment, size_t size) {
	return _aligned_malloc(size, alignment);
	
}
#else
inline void _aligned_free(void* _block)
{
	free(_block);
}
#endif