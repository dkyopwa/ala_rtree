﻿#include <stdio.h>
//#include <malloc.h>
#include <math.h>
//#include <stdlib.h>
#include <stdbool.h>
#ifndef _WIN
	#include <sys/time.h>
	#include <unistd.h>
	#include <pthread.h>
	#define __int64 long long
#else
	#include "time.h"
	#include <conio.h>
	#include <windows.h>  
	#include <process.h>    /* _beginthread, _endthread */  
#endif
#include "unimem.h"
#include "first.h"

struct test_data {
	struct node *nd;
	coord *xx;
	coord *yy;
	unsigned count;
	coord radius;
	indexer *idxs;
};

void find_test1(struct node *nd, const coord *xx, const coord *yy, const unsigned count, const coord radius, indexer *idxs);
void find_test2(struct node *nd, const coord *xx, const coord *yy, const unsigned count, const coord radius, indexer *idxs);
#ifndef _WIN
void* find_thread(void *params);
#else
void find_thread(void *params);
#endif
/// generate data for test
struct leaf* generate(unsigned *count, unsigned **offsets_leafs, unsigned *count_shapes);
struct leaf* automatic_generate(/*out*/unsigned *count, /*out*/unsigned **offsets_leafs, /*out*/unsigned *count_shapes, /*in*/unsigned max_count_shapes = 10000000, /*in*/unsigned max_points_in_shapes = 10);

void try_find(struct node *nd, struct leaf* lll, unsigned count_of_leafs);
void try_find2(struct node *nd, struct leaf* lll, unsigned count_of_leafs);
void try_find3(struct node *nd, struct leaf* lll, unsigned count_of_leafs);
void try_find4(struct node *nd, struct leaf* lll, unsigned count_of_leafs);
#ifdef USE_CUDA
#if defined(CALC_CIRCLE) || defined(CALC_POINT)
/* searchin items in selected rectangle on cuda device */
extern "C"
indexer* cuda_search_rect2(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, bool intersection, /*out*/indexer *count_items);
/* searching the nearest item to point in selected radius */
extern "C"
indexer* cuda_search_nearest_item2(/*in*//*struct node *nd, *//*in*/coord x, /*in*/coord y, /*in*/coord radius, bool intersection, /*out*/coord *dist);
#endif // CALC_POINT
void try_find5_cuda(struct node *nd, struct leaf* lll, unsigned count_of_leafs);
void try_find6_cuda(struct node *nd, struct leaf* lll, unsigned count_of_leafs);
#endif // USE_CUDA

// ---------------------------------------- TEST -------------------------------------------

int main()
{
	lprintf("Hello");

	unsigned count_of_leafs = 0; // old: must to devide 3 without
	alignas(16) unsigned *offsets_leafs = NULL;
	unsigned count_shapes = 0;
	struct leaf* lll = automatic_generate(&count_of_leafs, &offsets_leafs, &count_shapes);
	//struct leaf* lll = generate(&count_of_leafs, &offsets_leafs, &count_shapes);

	char ch[1024];
#ifdef _WIN
	sprintf_s(ch, 1024, "Shapes = %u, leafs = %u\n", count_shapes, count_of_leafs);
#else
	snprintf(ch, 1024, "Shapes = %u, leafs = %u\n", count_shapes, count_of_leafs);
#endif // WIN
	lprintf(ch);

	struct node* nd = create_rtree(lll, count_of_leafs, offsets_leafs, count_shapes);

	if (nd) {
		// testing
		lprintf("start");
		////try_find(nd, lll, count_of_leafs);
		//try_find2(nd, lll, count_of_leafs);
		////try_find3(nd, lll, count_of_leafs);
		//try_find4(nd, lll, count_of_leafs);
#ifdef USE_CUDA
		//try_find5_cuda(nd, lll, count_of_leafs);
		try_find6_cuda(nd, lll, count_of_leafs);
#endif //USE_CUDA
		lprintf("end");
#ifdef CALC_POINT
		//lprintf("=== POINT ===");
#endif // CALC_POINT
		}

	if (lll)
		_aligned_free(lll);
	_aligned_free(offsets_leafs);
	del_root();

	lprintf("Done");

#ifdef _WIN
	_getch();
#endif

	return 0;
}

struct leaf* automatic_generate(/*out*/unsigned *count, /*out*/unsigned **offsets_leafs, /*out*/unsigned *count_shapes, /*in*/unsigned max_count_shapes, /*in*/unsigned max_points_in_shapes)
{
	unsigned b = (unsigned)sqrt(max_count_shapes / 2.0); //707; // 707 2236
	unsigned a = b * 2;
	coord len = 360.0 / b / 4.0;
	*count_shapes = a * b;
	unsigned total_count = 0;
	unsigned tcount = 0;
	indexer big_number = 1;

	struct leaf *ret_leafs = (struct leaf*)aligned_alloc(16, sizeof(struct leaf) * a * b * (max_points_in_shapes + 1));
	*offsets_leafs = (unsigned*)aligned_alloc(16, sizeof(unsigned) * a * b);

	//struct data *dt = (struct data*)malloc(sizeof(struct data) * count);
	//struct point tpts;
	coord x, y, xc, yc;
	coord delta_x = 359.9 / a, delta_y = 179.9 / b;
	//unsigned idx;

	for (unsigned i = 0; i < a; ++i) {
		xc = -179.95 + delta_x * i;
		for (unsigned j = 0; j < b; ++j) {
			//idx = i + j * a;
			tcount = rand() % (max_points_in_shapes - 1) + 1;

			//dt[idx].count = tcount;
			//dt[idx].pts = (struct point*)malloc(sizeof(struct point) * tcount);

			for (unsigned k = 0; k < tcount; ++k) {
				yc = -89.95 + delta_y * j;

				coord alpha = (360.0 / tcount) * k * 2 * PI / 360.0;
				coord tlen = rand() % 2 + 1;
				x = xc + tlen * cos(alpha);
				y = yc + tlen * sin(alpha);

				ret_leafs[total_count].x = x;
				ret_leafs[total_count].y = y;
				ret_leafs[total_count].number = big_number;

				total_count++;
			}
			//if (tcount < 2)
			//	printf("BIG = %u, tcount = %u\n", big_number, tcount);
			(*offsets_leafs)[j + i * b] = tcount;
			big_number++;
		}
	}
	*count = total_count;

	ret_leafs = (struct leaf*)_aligned_realloc(ret_leafs, total_count * sizeof(struct leaf), 16);
	return ret_leafs;
}

/// generate test data
struct leaf* generate(unsigned *count, unsigned **offsets_leafs, unsigned *count_shapes)
{
	FILE *f1;
#ifndef _WIN
	f1 = fopen("/media/vovan/OS/projects/tmp/1/4.bin", "rb");
#else
	errno_t t = fopen_s(&f1, "c:/projects/tmp/1/4.bin", "rb");
#endif

	unsigned count1 = 0;
	unsigned count_leafs = 0;
	unsigned cnt = 0;
	size_t st = 0;
	unsigned ii = 0;
	double d[2];

	st = fread(&count1, sizeof(unsigned), 1, f1);
	st = fread(&count_leafs, sizeof(unsigned), 1, f1);
	printf("%u, leafs = %u\n", count1, count_leafs);
	*count = count_leafs;
	*offsets_leafs = (unsigned*)aligned_alloc(16, sizeof(unsigned) * count1);
	*count_shapes = count1;

	// prepare for ramdom number
	alignas(16) unsigned *numbers = (unsigned*)aligned_alloc(16, sizeof(unsigned) * count1);
	unsigned t2 = 0;
	unsigned total_number = 0;
	for (unsigned i = 0; i < count1; ++i) {
		t2 = (rand() * rand()) % count1;
		for (unsigned j = 0; j < total_number; ++j) {
			if (numbers[j] == t2) {
				t2++;
				if (t2 >= count1)
					t2 = 0;
				j = 0;
				continue;
			}
		}
		numbers[i] = t2;
	}

	alignas(16) struct leaf* res = NULL;
	res = (struct leaf*)aligned_alloc(16, sizeof(struct leaf) * count_leafs);
	if (!res)
		return NULL;

	// struct data *dts = (struct data*)malloc(sizeof(struct data) * count1);
	unsigned t1 = 0;
	for (unsigned i = 0; i < count1; ++i) {
		st = fread(&cnt, sizeof(unsigned), 1, f1);
		(*offsets_leafs)[i] = cnt;
		for (unsigned j = 0; j < cnt; ++j) {
			res[ii].number = i;
			st = fread(d, sizeof(double), 2, f1);
			res[ii].x = (coord)d[0];
			res[ii].y = (coord)d[1];
			++ii;
		}
		t1 += cnt;
		//dts[i].pts = (struct point*)malloc(sizeof(struct point) * cnt);
		//st = fread(dts[i].pts, sizeof(double), cnt * 2, f1);
	}
	printf("offsets sum = %u\n", t1);

	fclose(f1);
	_aligned_free(numbers);

	return res;
	// old variant
	//	struct leaf* res = NULL;
	/*	res = (struct leaf*)malloc(sizeof(struct leaf) * count);
	if (!res)
	return NULL;

	srand((unsigned int)time(NULL));
	for (unsigned i = 0; i < count; i += 3) {
	res[i].x = rand() % 9999 / 100.0;
	res[i].y = rand() % 9999 / 100.0;
	res[i].number = (unsigned)ceil(i / 3.0);
	coord alpha = ((coord)(rand() % 48) + 1.0) / 100.0 * PI;
	res[i + 1].x = res[i].x + (coord)((rand() % 5) / 100.0) / sin(alpha);
	res[i + 1].y = res[i].y + (coord)((rand() % 5) / 100.0) / cos(alpha);
	res[i + 1].number = (unsigned)ceil(i / 3.0);
	alpha = ((coord)(rand() % 48) + 1.0) / 100.0 * PI;
	res[i + 2].x = res[i].x + (coord)((rand() % 5) / 100.0) / sin(alpha);
	res[i + 2].y = res[i].y + (coord)((rand() % 5) / 100.0) / cos(alpha);
	res[i + 2].number = (unsigned)ceil(i / 3.0);

	char tch[1024] = {0};
	/////////		sprintf(tch, "%u: %f:%f, %f:%f, %f:%f", (unsigned)ceil(i / 3.0), res[i].x, res[i].y, res[i + 1].x, res[i + 1].y, res[i + 2].x, res[i + 2].y);
	lprintf(tch);
	}
	return res;
	*/
}

void try_find(struct node *nd, struct leaf* lll, unsigned count_of_leafs)
{
	/*	char **colors = (char**)malloc(sizeof(char*) * 125);
	for (unsigned i = 0; i < 125; ++i) {
	colors[i] = (char*)malloc(sizeof(char) * 32);
	}
	short t2 = 0;
	for (char i = 0; i < 5; ++i) {
	for (char j = 0; j < 5; ++j) {
	for (char k = 0; k < 5; ++k) {
	#ifndef _WIN
	snprintf(colors[t2], 32, "rgb(%u,%u,%u)", i * 51, j * 51, k * 51);
	#else
	sprintf_s(colors[t2], 32, "rgb(%u,%u,%u)", i * 51, j * 51, k * 51);
	#endif
	t2++;
	}
	}
	}
	*/
	const unsigned count = 1000000;
	/*const coord radius = (nd->x2 - nd->x1) < (nd->y2 - nd->y1) ? (nd->x2 - nd->x1)/100.0 : (nd->y2 - nd->y1)/100.0;
	char ch1[1024];
	sprintf_s(ch1, 1024, "RADIUS = %f", radius);
	lprintf(ch1);
	_getch();*/
	const coord radius = (coord)0.3;

	alignas(16) indexer *idxs1 = (indexer*)aligned_alloc(16, sizeof(indexer) * count);
	alignas(16) indexer *idxs2 = (indexer*)aligned_alloc(16, sizeof(indexer) * count);
	alignas(16) coord *xx = (coord*)aligned_alloc(16, sizeof(coord) * count);
	alignas(16) coord *yy = (coord*)aligned_alloc(16, sizeof(coord) * count);
	for (unsigned i = 0; i < count; ++i) {
		xx[i] = (coord)((rand() % ((int)((nd->x2 - nd->x1) * 10.0))) / 10.0 + nd->x1); // - (nd->x2 - nd->x1) / 2.0;
		yy[i] = (coord)((rand() % ((int)((nd->y2 - nd->y1) * 10.0))) / 10.0 + nd->y1); //- (nd->y2 - nd->y1) / 2.0;
	}

	lprintf("start 1");
	find_test1(nd, xx, yy, count, radius, idxs1);
	lprintf("end 1");
	lprintf("start 2");
	find_test2(nd, xx, yy, count, radius, idxs2);
	lprintf("end 2");
	/*	indexer idx;
	char ch[64];
	lprintf("start2");
	for (unsigned i = 0; i < count; i++) {
	idx = search_point(m_nodes, xx[i], yy[i], radius);
	idxs[i] = idx;
	if (idx != (indexer)-1) {
	//sprintf_s(ch, 64, "Find = %u", idx);
	//lprintf(ch);
	} else {
	#ifndef _WIN
	snprintf(ch, 64, "Error %u (x = %u, y = %u)", i, (unsigned)xx[i], (unsigned)yy[i]);
	#else
	sprintf_s(ch, 64, "Error %u (x = %u, y = %u)", i, (unsigned)xx[i], (unsigned)yy[i]);
	#endif
	lprintf(ch);
	}
	}
	*/

	// compare results
	for (indexer i = 0; i < count; ++i) {
		if (idxs1[i] != idxs2[i]) {
			printf("Error %u: %u vs %u\n", i, idxs1[i], idxs2[i]);
		}
	}

	/*	FILE *f2;
	char name[256];
	sprintf_s(name, 256, "c:/projects/tmp/1/%s", "test_res.svg");
	errno_t t = fopen_s(&f2, name, "w");
	unsigned t1 = 0;
	short num_color = 1;

	fprintf(f2, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" height=\"10000px\"  width=\"10000px\">\n"); //  height=\"400px\"  width=\"400px\"
	for (unsigned k = 0; k < count; ++k) {
	unsigned i;
	if (idxs[k] == (indexer)-1)
	continue;
	for (i = 0; i < count_of_leafs; ++i) {
	if (lll[i].number == idxs[k])
	break;
	}
	if (i >= count_of_leafs)
	continue;
	fprintf(f2, "\t<polygon points=\"");
	while (lll[i].number == idxs[k]) {
	fprintf(f2, "%u,%u ", (unsigned)lll[i].x, (unsigned)lll[i].y);
	i++;
	}
	//fprintf(f2, "%u,%u ", (unsigned)xx[k], (unsigned)yy[k]);
	fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(0, 0, 0)\" fill=\"%s\"/>\n", colors[num_color]); // rgb(0, 0, 0) rgb(150,150,255)
	fprintf(f2, "\t<circle cx=\"%u\" cy=\"%u\" r=\"3\" stroke-width=\"1\" stroke=\"rgb(50, 50, 50)\" fill=\"%s\"/>\n", (unsigned)xx[k], (unsigned)yy[k], colors[num_color]);
	num_color++;
	if (num_color > 124)
	num_color = 1;
	}

	fprintf(f2, "</svg>");
	fclose(f2);
	*/
	_aligned_free(xx);
	_aligned_free(yy);
	_aligned_free(idxs1);
	_aligned_free(idxs2);
	/*	for (unsigned i = 0; i < 125; ++i) {
	free(colors[i]);
	}
	free(colors);
	*/
}

void find_test1(struct node *nd, const coord *xx, const coord *yy, const unsigned count, const coord radius, indexer *idxs)
{
	indexer idx;
	char ch[64];
	//lprintf("start2");
	for (unsigned i = 0; i < count; i++) {
		idx = search_point(nd, xx[i], yy[i], radius);
		idxs[i] = idx;
		if (idx != (indexer)-1) {
			//sprintf_s(ch, 64, "Find = %u", idx);
			//lprintf(ch);
		}
		else {
#ifndef _WIN
			snprintf(ch, 64, "Error %u (x = %u, y = %u)", i, (unsigned)xx[i], (unsigned)yy[i]);
#else
			sprintf_s(ch, 64, "Error %u (x = %u, y = %u)", i, (unsigned)xx[i], (unsigned)yy[i]);
#endif
			lprintf(ch);
		}
	}
}

void find_test2(struct node *nd, const coord *xx, const coord *yy, const unsigned count, const coord radius, indexer *idxs)
{
	unsigned cpus = 1;
#ifndef _WIN
	cpus = sysconf(_SC_NPROCESSORS_ONLN);
#else
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	cpus = sysinfo.dwNumberOfProcessors;
#endif
	//cpus = 1; /*!!!!!!!*/

#ifndef _WIN
	// handler thread
	pthread_t *ptr1 = (pthread_t*)malloc(sizeof(pthread_t) * cpus);
	//unsigned tun[32];
	struct test_data *data = (struct test_data*)malloc(sizeof(struct test_data) * cpus);
	//unsigned *idx = (unsigned*)malloc(sizeof(unsigned) * (count_cpus + 1));

	// prepare for separate
	unsigned offset = (unsigned)ceil(count / cpus);
	//idx[0] = 0;
	//idx[count_cpus] = count_leafs;

	for (unsigned i = 0; i < cpus; ++i) {
		data[i].count = offset;
		if (i == cpus - 1) {
			data[i].count = count % offset;
			if (!data[i].count)
				data[i].count = offset;
		}

		data[i].idxs = idxs + offset * i;
		data[i].nd = nd;
		data[i].radius = radius;
		data[i].xx = (coord*)xx + offset * i;
		data[i].yy = (coord*)yy + offset * i;

		//tun[i] = i;
		pthread_create(&(ptr1[i]), NULL, find_thread, &(data[i]));
	}
	for (unsigned i = 0; i < cpus; ++i) {
		//WaitForSingleObject((HANDLE)ptr1[i], INFINITE);
		pthread_join(ptr1[i], NULL);
	}

	//free(idx);
	free(data);
	free(ptr1);

#else
	// handler thread
	alignas(16) uintptr_t *ptr1 = (uintptr_t*)aligned_alloc(16, sizeof(uintptr_t) * cpus);
	//unsigned tun[32];
	alignas(16) struct test_data *data = (struct test_data*)aligned_alloc(16, sizeof(struct test_data) * cpus);
	//unsigned *idx = (unsigned*)malloc(sizeof(unsigned) * (cpus + 1));

	// prepare for separate
	unsigned offset = (unsigned)ceil((double)count / cpus);
	//idx[0] = 0;
	//idx[count_cpus] = count_leafs;
	
	for (unsigned i = 0; i < cpus; ++i) {
		data[i].count = offset;
		if (i == cpus - 1) {
			data[i].count = count % offset;
			if (!data[i].count)
				data[i].count = offset;
		}

		data[i].idxs = idxs + offset * i;
		data[i].nd = nd;
		data[i].radius = radius;
		data[i].xx = (coord*)xx + offset * i;
		data[i].yy = (coord*)yy + offset * i;

		//tun[i] = i;
		ptr1[i] = _beginthread(find_thread, 0, &(data[i]));
	}
	for (unsigned i = 0; i < cpus; ++i) {
		WaitForSingleObject((HANDLE)ptr1[i], INFINITE);
	}

	//free(idx);
	_aligned_free(data);
	_aligned_free(ptr1);
#endif // _WIN
}

#ifndef _WIN
void* find_thread(void *params)
#else
void find_thread(void *params)
#endif
{
	struct test_data *data = (struct test_data*)params;
	find_test1(data->nd, data->xx, data->yy, data->count, data->radius, data->idxs);
#ifndef _WIN
	return NULL;
#else
	return;
#endif
}

void try_find2(struct node *nd, struct leaf* lll, unsigned count_of_leafs)
{
	indexer *idxs1 = NULL, *idxs2 = NULL, *idxs3 = NULL, *idxs4 = NULL;
	indexer count1 = 0, count2 = 0, count3 = 0, count4 = 0;
	lprintf("start3");
#ifdef _WIN
	__int64 t1 = __rdtsc();
#else
	clock_t t1 = clock();
#endif //_win
	idxs1 = search_in_rect(nd, 50, 50, 55, 55, &count1);
#ifdef _WIN
	__int64 t2 = __rdtsc();
#else
	clock_t t2 = clock();
#endif //_WIN
	idxs2 = search_rect2(nd, 50, 50, 55, 55, false, &count2);
	//idxs4 = search_rect2(nd, 50, 50, 55, 55, true, &count4, ret_callback2_circle(1));
#ifdef _WIN
	__int64 t3 = __rdtsc();
#else
	clock_t t3 = clock();
#endif //_WIN
#ifdef CALC_CIRCLE
	idxs3 = search_circle2(nd, 52.5, 52.5, 2.5, true, &count3);
#elif defined(CALC_POINT)
	coord dist = 0;
	idxs3 = search_nearest_item2(nd, 52.5, 52.5, 2.5, false, &count3, &dist);
#endif // CALC_CIRCLE
	char ch[1024];
#ifdef _WIN
	__int64 t4 = __rdtsc();
	sprintf_s(ch, 1024, "end3 %lld vs %lld => %f (c1 = %u, c2 = %u, c4 = %u)\nend4 time = %lld (count = %d) speed = %f", t2 - t1, t3 - t2, (t3 - t2) * 100.0 / (t2 - t1), count1, count2, count4, t4 - t3, count3, (t4 - t3) * 100.0 / (t3 - t2));
#else
	clock_t t4 = clock();
	snprintf(ch, 1024, "end3 %ld vs %ld => %f (c1 = %u, c2 = %u, c4 = %u)\nend4 time = %ld (count = %d) speed = %f", t2 - t1, t3 - t2, (t3 - t2) * 100.0 / (t2 - t1), count1, count2, count4, t4 - t3, count3, (t4 - t3) * 100.0 / (t3 - t2));
#endif //_WIN
	lprintf(ch);
#ifdef CALC_POINT
#ifdef _WIN
	sprintf_s(ch, 1024, "Near item dist = %f, idx = %u", dist, idxs3[0]);
#else
	snprintf(ch, 1024, "Near item dist = %f, idx = %u", dist, idxs3[0]);
#endif //_WIN
	lprintf(ch);
#endif
	if (idxs1)
		_aligned_free(idxs1);
	if (idxs2)
		_aligned_free(idxs2);
	if (idxs3)
		_aligned_free(idxs3);
	if (idxs4)
		_aligned_free(idxs4);
/*	if (idxs) {
		printf("count = %u\n", count);

		FILE *f2;
		char name[256];
		sprintf_s(name, 256, "c:/projects/tmp/1/%s", "test2_res.svg");
		errno_t t = fopen_s(&f2, name, "w");
		//unsigned t1 = 0;
		//short num_color = 1;

		fprintf(f2, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" height=\"10000px\"  width=\"10000px\">\n"); //  height=\"400px\"  width=\"400px\"
		fprintf(f2, "\t<polygon points=\"100,100 1012,100 1012,1000 100,1000\" stroke-width=\"1\" stroke=\"rgb(50, 50, 150)\" fill=\"none\"/>\n");
		*/
		/*for (unsigned k = 0; k < count; ++k) {
			unsigned i;
			if (idxs[k] == (indexer)-1)
				continue;
			for (i = 0; i < count_of_leafs; ++i) {
				if (lll[i].number == idxs[k])
					break;
			}
			if (i >= count_of_leafs)
				continue;
			fprintf(f2, "\t<polygon points=\"");
			while (lll[i].number == idxs[k]) {
				fprintf(f2, "%u,%u ", (unsigned)lll[i].x, (unsigned)lll[i].y);
				i++;
			}
			//fprintf(f2, "%u,%u ", (unsigned)xx[k], (unsigned)yy[k]);
			fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(0, 0, 0)\" fill=\"%s\"/>\n", "rgb(150, 150, 255)"); // rgb(0, 0, 0) rgb(150,150,255)
			//fprintf(f2, "\t<circle cx=\"%u\" cy=\"%u\" r=\"3\" stroke-width=\"1\" stroke=\"rgb(50, 50, 50)\" fill=\"%s\"/>\n", (unsigned)xx[k], (unsigned)yy[k], colors[num_color]);
		} */

/*		for (indexer i = 0; i < count_of_leafs; ++i) {
			bool in_rect = false;
			indexer tmp = lll[i].number;
			for (indexer k = 0; k < count; ++k) {
				if (idxs[k] == tmp) {
					in_rect = true;
					break;
				}
			}
			fprintf(f2, "\t<polygon points=\"");
			while (lll[i].number == tmp) {
				fprintf(f2, "%u,%u ", (unsigned)lll[i].x, (unsigned)lll[i].y);
				i++;
			}
			if (in_rect) {
				fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(0, 0, 0)\" fill=\"%s\"/>\n", "rgb(150, 150, 255)");
			}
			else {
				fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(0, 0, 0)\" fill=\"%s\"/>\n", "rgb(255, 150, 150)");
			}
		}

		fprintf(f2, "</svg>");
		fclose(f2);
		
	}*/

//	if (idxs)
//		_aligned_free(idxs);
}

void try_find3(struct node *nd, struct leaf* lll, unsigned count_of_leafs)
{
	alignas(16) indexer *idxs = NULL;
	indexer count = 0;
	lprintf("start4");
	idxs = search_in_circles(nd, 1000, 1000, 500.0, &count);
	lprintf("end4");
	if (idxs) {
		printf("count = %u\n", count);

/*		FILE *f2;
		char name[256];
		sprintf_s(name, 256, "c:/projects/tmp/1/%s", "test3_res.svg");
		errno_t t = fopen_s(&f2, name, "w");
		//unsigned t1 = 0;
		//short num_color = 1;

		fprintf(f2, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" height=\"10000px\"  width=\"10000px\">\n"); //  height=\"400px\"  width=\"400px\"
		fprintf(f2, "\t<circle cx=\"%u\" cy=\"%u\" r=\"500\" stroke-width=\"1\" stroke=\"rgb(50, 255, 100)\" fill=\"none\"/>\n", 1000, 1000);
		for (unsigned k = 0; k < count; ++k) {
		unsigned i;
		if (idxs[k] == (indexer)-1)
		continue;
		for (i = 0; i < count_of_leafs; ++i) {
		if (lll[i].number == idxs[k])
		break;
		}
		if (i >= count_of_leafs)
		continue;
		fprintf(f2, "\t<polygon points=\"");
		while (lll[i].number == idxs[k]) {
		fprintf(f2, "%u,%u ", (unsigned)lll[i].x, (unsigned)lll[i].y);
		i++;
		}
		//fprintf(f2, "%u,%u ", (unsigned)xx[k], (unsigned)yy[k]);
		fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(0, 0, 0)\" fill=\"%s\"/>\n", "rgb(150, 150, 255)"); // rgb(0, 0, 0) rgb(150,150,255)
		//fprintf(f2, "\t<circle cx=\"%u\" cy=\"%u\" r=\"3\" stroke-width=\"1\" stroke=\"rgb(50, 50, 50)\" fill=\"%s\"/>\n", (unsigned)xx[k], (unsigned)yy[k], colors[num_color]);
		}

		fprintf(f2, "</svg>");
		fclose(f2);
		*/
	}

	if (idxs)
		_aligned_free(idxs);
}

void try_find4(struct node *nd, struct leaf* lll, unsigned count_of_leafs)
{
	/*float *fl_x = (float*)malloc(sizeof(float) * 3600 * 1800);
	float *fl_y = (float*)malloc(sizeof(float) * 3600 * 1800);
	for (float xi = -179.9; xi < 180.0; xi += 0.1) {
		for (float yi = -89.9; yi < 90.0; yi += 0.1) {
			fl_x[xi * 3600 + yi] = xi;
			fl_y[xi * 3600 + yi] = yi;
		}
	}*/
	indexer count_items;
	coord dist, dd, mind, maxd;
	unsigned i = 0;
	clock_t t1 = clock();
	for (float yi = -89.9; yi < 90.0; yi += 0.1, ++i) {
		clock_t t11 = clock();
		dd = 0.0;
		mind = 9999999999.9;
		maxd = 0.0;
		for (float xi = -180.9; xi < 180.5; xi += 0.1) {
			indexer *idxs = search_nearest_item2(nd, xi, yi, 1.0, true, &count_items, &dist);
			dd += dist;
			if (count_items > 1)
				printf("Error %f, %f\n", xi, yi);
			if (dist > maxd)
				maxd = dist;
			if (dist < mind)
				mind = dist;
			_aligned_free(idxs);
		}
		clock_t t22 = clock();
		printf("Line %u: time = %u ms, min dist = %f, max dist = %f, avg dist = %f\n", i, t22 - t11, mind, maxd, dd / 3598.0);
	}
	clock_t t2 = clock();


	/* free(fl_y);
	free(fl_x); */
}

void try_find5_cuda(struct node *nd, struct leaf* lll, unsigned count_of_leafs)
{
	int count_iter = 100;

	indexer count_items1 = 0, count_items2 = 0;
	//coord dist;
	clock_t t1 = clock();
	indexer *idxs1 = NULL;
	for (int i = 0; i < count_iter; ++i) {
		idxs1 = search_rect2(nd, -50, 1, 1, 2.0, false, &count_items1); //-180, -87, 180, -86.6
		if (i != count_iter - 1)
			_aligned_free(idxs1);
	}
	clock_t t2 = clock();
	printf("Time to calculate on CPU: %d, iter = %i\n", t2 - t1, count_iter);
#ifdef USE_CUDA
	t1 = clock();
	indexer *idxs2 = NULL;
	for (int i = 0; i < count_iter; ++i) {
		idxs2 = cuda_search_rect2(nd, -50, 1, 1, 2.0, false, &count_items2);
		if (i != count_iter - 1)
			_aligned_free(idxs2);
	}
	t2 = clock();
	printf("Time to calculate on GPU: %d, iter = %i\n", t2 - t1, count_iter);
#endif
	if (count_items1 != count_items2 && count_items2 != 99999) {
		printf("Error in count items: %u vs %u\n", count_items1, count_items2);

		/*for (indexer i = 0; i < count_items2; ++i) {
			bool fl = false;
			for (indexer j = 0; j < count_items1; ++j) {
				if (idxs1[j] == idxs2[i]) {
					fl = true;
					break;
				}
			}
			if (!fl) {
				printf("%u\n", idxs2[i]);
			}
		}*/

/*		printf("RESULT 1\n");
		for (indexer i = 0; i < count_items1; ++i)
			printf("%u\n", idxs1[i]);
		unsigned start = 0, stop = 0;
		bool flag = false;
		for (indexer i = 0; i < count_items1; ++i) {
			start = 0; stop = 0;
			for (unsigned j = 0; j < count_of_leafs; ++j) {
				if (lll[j].number == idxs1[i]) {
					if (!start)
						start = j;
					stop = j + 1;
				}
				if (lll[j].number != idxs1[i] && start && stop)
					break;
			}
			//if (stop - 1 == start)
			//	printf("Error shape N %u\n", idxs1[i]);
			for (unsigned k = start; k < stop; ++k) {
				if (lll[k].x > -1 && lll[k].x < 1 && lll[k].y > 1 && lll[k].y < 1.5) {
					flag = true;
					break;
				}
			}
			if (!flag)
				for (unsigned k = start; k < stop; ++k)
					printf("IDX 1: %u (%u): x = %f, y = %f\n", idxs1[i], lll[k].number, lll[k].x, lll[k].y);
			flag = false;
		}
		*/
#ifdef USE_CUDA
	/*	printf("\n\n\nRESULT 2\n");
		for (indexer i = 0; i < count_items2; ++i)
			printf("%u\n", idxs2[i]);
		for (indexer i = 0; i < count_items2; ++i) {
			start = 0; stop = 0;
			for (unsigned j = 0; j < count_of_leafs; ++j) {
				if (lll[j].number == idxs2[i]) {
					if (!start)
						start = j;
					stop = j + 1;
				}
				if (lll[j].number != idxs2[i] && start && stop)
					break;
			}
			//if (stop - 1 == start)
			//	printf("Error shape N %u\n", idxs2[i]);
			for (unsigned k = start; k < stop; ++k) {
				if (lll[k].x > -1 && lll[k].x < 1 && lll[k].y > 1 && lll[k].y < 1.5) {
					flag = true;
					break;
				}
			}
			if (!flag)
				for (unsigned k = start; k < stop; ++k)
					printf("IDX 2: %u (%u): x = %f, y = %f\n", idxs2[i], lll[k].number, lll[k].x, lll[k].y);
			flag = false;
		}*/
#endif
	}
	else {
		printf("Result count = %u (%u)\n", count_items1, count_items2);
	}

#ifdef USE_CUDA
	_aligned_free(idxs2);
#endif
	_aligned_free(idxs1);
}

void try_find6_cuda(struct node *nd, struct leaf* lll, unsigned count_of_leafs)
{
	int count_iter = 1000;
	coord radius = 2.0;

	indexer *tidxs1 = (indexer*)malloc(sizeof(indexer) * count_iter);
	coord *tdist1 = (coord*)malloc(sizeof(coord) * count_iter);
	indexer *tidxs2 = (indexer*)malloc(sizeof(indexer) * count_iter);
	coord *tdist2 = (coord*)malloc(sizeof(coord) * count_iter);
	coord *x = (coord*)malloc(sizeof(coord) * count_iter);
	coord *y = (coord*)malloc(sizeof(coord) * count_iter);

	indexer *idxs1 = NULL, *idxs2 = NULL;
	coord dist1 = 0.0, dist2 = 0.0;
	indexer c1 = 0, c2 = 0;

	for (int i = 0; i < count_iter; ++i) {
		x[i] = (rand() % 3600) / 10.0 - 180.0; //91.9;
		y[i] = (rand() % 1800) / 10.0 - 90.0; //-55.4
	}

	clock_t t1 = clock();
	for (int i = 0; i < count_iter; ++i) {
		idxs1 = search_nearest_item2(nd, x[i], y[i], radius, true, &c1, &dist1);
		tidxs1[i] = idxs1[0];
		tdist1[i] = dist1;
		_aligned_free(idxs1);
		//printf("CPU result = %u, dist = %f\n", tidxs1[i], tdist1[i]);
	}
	clock_t t2 = clock();
	printf("Time to search the nearest item on CPU = %d, count = %i\n", t2 - t1, count_iter);

#ifdef USE_CUDA
	clock_t t3 = clock();
	for (int i = 0; i < count_iter; ++i) {
		idxs2 = cuda_search_nearest_item2(/*nd, */x[i], y[i], radius, true, &dist2);
		tidxs2[i] = idxs2[0];
		tdist2[i] = dist2;
		_aligned_free(idxs2);
		//printf("GPU result = %u, dist = %f\n", tidxs2[i], tdist2[i]);
	}
	clock_t t4 = clock();
	printf("Time to search the nearest item on GPU = %d, count = %i\n", t4 - t3, count_iter);

	for (int i = 0; i < count_iter; ++i) {
		if (tidxs1[i] != tidxs2[i] && tdist1[i] != tdist2[i]) {
			printf("========== Error %i: idx = %u vs %u, dist = %e vs %e, coord = %f, %f\n", i, tidxs1[i], tidxs2[i], tdist1[i], tdist2[i], x[i], y[i]);
		}
	}
#endif

	/*if (idxs1)
		_aligned_free(idxs1);
	if (idxs2)
		_aligned_free(idxs2);
		*/
	free(x);
	free(y);
	free(tidxs1);
	free(tidxs2);
	free(tdist1);
	free(tdist2);
}