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
/// generate data for test
struct leaf* generate(unsigned *count, unsigned **offsets_leafs, unsigned *count_shapes);
void try_find(struct node *nd, struct leaf* lll, unsigned count_of_leafs);
void try_find2(struct node *nd, struct leaf* lll, unsigned count_of_leafs);
void try_find3(struct node *nd, struct leaf* lll, unsigned count_of_leafs);

// ---------------------------------------- TEST -------------------------------------------

int main()
{
	lprintf("Hello");

	unsigned count_of_leafs = 0; // old: must to devide 3 without
	unsigned *offsets_leafs = NULL;
	unsigned count_shapes = 0;
	struct leaf* lll = generate(&count_of_leafs, &offsets_leafs, &count_shapes);
	//lprintf("Done"); _getch();  return 0;

	struct node* nd = create_rtree(lll, count_of_leafs, offsets_leafs, count_shapes);

	if (nd) {
		// testing
		lprintf("start");
		try_find(nd, lll, count_of_leafs);
		try_find2(nd, lll, count_of_leafs);
		try_find3(nd, lll, count_of_leafs);
		lprintf("end");
	}

	if (lll)
		free(lll);
	free(offsets_leafs);
	del_root();

	lprintf("Done");

#ifdef _WIN
	_getch();
#endif

	return 0;
}

/// generate test data
struct leaf* generate(unsigned *count, unsigned **offsets_leafs, unsigned *count_shapes)
{
	FILE *f1;
#ifndef _WIN
	f1 = fopen("/media/vovan/OS/projects/tmp/1/7.bin", "rb");
#else
	errno_t t = fopen_s(&f1, "c:/projects/tmp/1/7.bin", "rb");
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
	*offsets_leafs = (unsigned*)malloc(sizeof(unsigned) * count1);
	*count_shapes = count1;

	// prepare for ramdom number
	unsigned *numbers = (unsigned*)malloc(sizeof(unsigned) * count1);
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

	struct leaf* res = NULL;
	res = (struct leaf*)malloc(sizeof(struct leaf) * count_leafs);
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
	free(numbers);

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
	const unsigned count = 50000;
	const coord radius = 100;

	indexer *idxs1 = (indexer*)malloc(sizeof(indexer) * count);
	indexer *idxs2 = (indexer*)malloc(sizeof(indexer) * count);
	coord *xx = (coord*)malloc(sizeof(coord) * count);
	coord *yy = (coord*)malloc(sizeof(coord) * count);
	for (unsigned i = 0; i < count; ++i) {
		xx[i] = (rand() % ((int)((nd->x2 - nd->x1) * 10))) / 10.0 + nd->x1; // - (nd->x2 - nd->x1) / 2.0;
		yy[i] = (rand() % ((int)((nd->y2 - nd->y1) * 10))) / 10.0 + nd->y1; //- (nd->y2 - nd->y1) / 2.0;
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
	uintptr_t *ptr1 = (uintptr_t*)malloc(sizeof(uintptr_t) * cpus);
	//unsigned tun[32];
	struct test_data *data = (struct test_data*)malloc(sizeof(struct test_data) * cpus);
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

void try_find2(struct node *nd, struct leaf* lll, unsigned count_of_leafs)
{
	indexer *idxs = NULL;
	indexer count = 0;
	lprintf("start3");
	idxs = search_in_rect(nd, 100, 100, 512, 500, &count);
	lprintf("end3");
	if (idxs) {
		printf("count = %u\n", count);

/*		FILE *f2;
		char name[256];
		sprintf_s(name, 256, "c:/projects/tmp/1/%s", "test2_res.svg");
		errno_t t = fopen_s(&f2, name, "w");
		//unsigned t1 = 0;
		//short num_color = 1;

		fprintf(f2, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" height=\"10000px\"  width=\"10000px\">\n"); //  height=\"400px\"  width=\"400px\"
		fprintf(f2, "\t<polygon points=\"100,100 512,100 512,500 100,500\" stroke-width=\"1\" stroke=\"rgb(50, 50, 150)\" fill=\"none\"/>\n");
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
		free(idxs);
}

void try_find3(struct node *nd, struct leaf* lll, unsigned count_of_leafs)
{
	indexer *idxs = NULL;
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
		free(idxs);
}
