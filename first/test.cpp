#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#ifndef _WIN
	#include <sys/time.h>
	#include <unistd.h>
	#include <pthread.h>
#else
	#include "time.h"
	#include <conio.h>
	#include <windows.h>  
	#include <process.h>    /* _beginthread, _endthread */  
#endif
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
// ---------------------------------------- TEST -------------------------------------------

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
	const coord radius = 100;

	indexer *idxs1 = (indexer*)malloc(sizeof(indexer) * count);
	indexer *idxs2 = (indexer*)malloc(sizeof(indexer) * count);
	coord *xx = (coord*)malloc(sizeof(coord) * count);
	coord *yy = (coord*)malloc(sizeof(coord) * count);
	for (unsigned i = 0; i < count; ++i) {
		xx[i] = rand() % (int)nd->x2;
		yy[i] = rand() % (int)nd->y2;
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
	free(xx);
	free(yy);
	free(idxs1);
	free(idxs2);
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
	//cpus = 1;

#ifndef _WIN
	// handler thread
	pthread_t *ptr1 = (pthread_t*)malloc(sizeof(pthread_t) * count_cpus);
	//unsigned tun[32];
	struct first_thr_st *ftst = (struct first_thr_st*)malloc(sizeof(struct first_thr_st) * count_cpus);
	unsigned *idx = (unsigned*)malloc(sizeof(unsigned) * (count_cpus + 1));

	// prepare for separate
	unsigned offset = (unsigned)ceil(count_leafs / count_cpus);
	idx[0] = 0;
	idx[count_cpus] = count_leafs;
	for (unsigned i = 1; i < count_cpus; ++i) {
		for (unsigned j = offset * i; j < offset * i + offset - 1; ++j) {
			if (leafs[j].number != leafs[j + 1].number) {
				idx[i] = j + 1;
				break;
			}
		}
	}

	for (unsigned i = 0; i < count_cpus; ++i) {
		ftst[i].leafs_ = leafs;
		ftst[i].node_ = &(nd[i]);
		ftst[i].offsets_leafs_ = &(offsets_leafs[leafs[idx[i]].number]);
		ftst[i].count_leaf = idx[i + 1] - idx[i];
		ftst[i].start_pos_leafs = idx[i];

		//tun[i] = i;
		//ptr1[i] = _beginthread(first_thread_v2, 0, &(ftst[i]));
		pthread_create(&(ptr1[i]), NULL, first_thread_v2, &(ftst[i]));
	}
	for (unsigned i = 0; i < count_cpus; ++i) {
		//WaitForSingleObject((HANDLE)ptr1[i], INFINITE);
		pthread_join(ptr1[i], NULL);
	}

	free(idx);
	free(ftst);
	free(ptr1);

#else
	// handler thread
	uintptr_t *ptr1 = (uintptr_t*)malloc(sizeof(uintptr_t) * cpus);
	//unsigned tun[32];
	struct test_data *data = (struct test_data*)malloc(sizeof(struct test_data) * cpus);
	//unsigned *idx = (unsigned*)malloc(sizeof(unsigned) * (cpus + 1));

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
		ptr1[i] = _beginthread(find_thread, 0, &(data[i]));
	}
	for (unsigned i = 0; i < cpus; ++i) {
		WaitForSingleObject((HANDLE)ptr1[i], INFINITE);
	}

	//free(idx);
	free(data);
	free(ptr1);
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
