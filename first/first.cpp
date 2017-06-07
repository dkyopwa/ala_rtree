#include "pch.h"
#include <stdio.h>
#include <malloc.h>
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
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "sort.h"
#include "log.h"
#include "first.h"

//#define PRINT_SVG

/*#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h> // for close
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <errno.h>
#include <time.h>
#include <poll.h>
*/

// globale variables
struct node* m_nodes = NULL;
//struct center_node_st* m_branch_center = NULL;

/// main function
int main()
{
	lprintf("Hello");

	unsigned count_of_leafs = 0; // old: must to devide 3 without
	unsigned *offsets_leafs = NULL;
	unsigned count_shapes = 0;
	struct leaf* lll = generate(&count_of_leafs, &offsets_leafs, &count_shapes);
	//lprintf("Done"); _getch();  return 0;

	// threads
	unsigned cpus = 8;
#ifndef _WIN
	cpus = sysconf(_SC_NPROCESSORS_ONLN);
	char tch[1024];
	snprintf(tch, 1024, "CPUS = %u", cpus);
#else
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	cpus = sysinfo.dwNumberOfProcessors;
	char tch[1024];
	sprintf_s(tch, 1024, "CPUS = %u", cpus);
#endif
	lprintf(tch);

	struct node *nds = (struct node*)malloc(sizeof(struct node) * cpus);
	create_first_thread(nds, lll, offsets_leafs, cpus, count_of_leafs);
	lprintf("Print file");
	// print file in svg format only for 6.bin
	//print_file_svg(nds, cpus, 0, "7_.svg");

	// copy memory for united all leaf in one node
	// think about exclude first step and add all leafs in one node and calc only boundary
	unsigned total_mem = 0;
	for (unsigned i = 0; i < cpus; ++i) {
		total_mem += ((struct branch*)(nds[i].child_node[0]))->count_leafs;
	}
	total_mem = (unsigned)(ceil(total_mem / 16.0) * 16);
	//((struct branch*)(nds[0].child_node[0]))->leafs = (struct leaf*)realloc(((struct branch*)(nds[0].child_node[0]))->leafs, sizeof(struct leaf) * total_mem);
	init_root(&(lll[0]));
	struct branch *br1 = (struct branch*)(m_nodes->child_node[0]);
	br1->leafs = (struct leaf*)malloc(sizeof(struct leaf) * total_mem);
	br1->merge_next_leaf = (bool*)malloc(sizeof(bool) * total_mem);
	for (unsigned i = 0; i < cpus; ++i) {
		// copy leafs
		memcpy(br1->leafs + br1->count_leafs, ((struct branch*)(nds[i].child_node[0]))->leafs, ((struct branch*)(nds[i].child_node[0]))->count_leafs * sizeof(struct leaf));
		memcpy(br1->merge_next_leaf + br1->count_leafs, ((struct branch*)(nds[i].child_node[0]))->merge_next_leaf, ((struct branch*)(nds[i].child_node[0]))->count_leafs * sizeof(bool));
		br1->count_leafs += ((struct branch*)(nds[i].child_node[0]))->count_leafs;
		// find max boundary
		if (m_nodes->x1 > nds[i].x1)
			m_nodes->x1 = nds[i].x1;
		if (m_nodes->y1 > nds[i].y1)
			m_nodes->y1 = nds[i].y1;
		if (m_nodes->x2 < nds[i].x2)
			m_nodes->x2 = nds[i].x2;
		if (m_nodes->y2 < nds[i].y2)
			m_nodes->y2 = nds[i].y2;
		
		// free memory
		free(((struct branch*)(nds[i].child_node[0]))->merge_next_leaf);
		free(((struct branch*)(nds[i].child_node[0]))->leafs);
		free(((struct branch*)(nds[i].child_node[0])));
		free(nds[i].center_child_node);
	}

	/*FILE *f3 = NULL;
	fopen_s(&f3, "c:/projects/tmp/1/tmp.txt", "w");

	for (unsigned i = 0; i < br1->count_leafs; ++i) {
		fprintf(f3, "%u:\tnum=%u\t%u\n", i, br1->leafs[i].number, (unsigned)(br1->merge_next_leaf[i]));
	}
	fclose(f3);

	_getch();
	return 0;
	*/

	// find centers
	if (!find_centers(m_nodes, count_shapes)) {
		lprintf("Centers error");
		print_file_svg(m_nodes, 1, 0, "7_0.svg");
	} else {
		print_file_svg(m_nodes, 1, count_shapes, "7_0.svg");
		separate(m_nodes);

		// calculate boundary of branches and lengths
		for (indexer i = 0; i < m_nodes->count_child_nodes; ++i) {
			if (!m_nodes->is_last_node) {
				lprintf("TO DO with first step separate branches");
				break;
			} 

			struct branch *tbr = (struct branch*)(m_nodes->child_node)[i];
			tbr->length = (coord*)malloc(sizeof(coord) * tbr->count_leafs);

			tbr->x_min = tbr->leafs[0].x;
			tbr->x_max = tbr->leafs[0].x;
			tbr->y_min = tbr->leafs[0].y;
			tbr->y_max = tbr->leafs[0].y;
			for (indexer j = 1; j < tbr->count_leafs; ++j) {
				// boundary
				if (tbr->leafs[j].x < tbr->x_min)
					tbr->x_min = tbr->leafs[j].x;
				else if (tbr->leafs[j].x > tbr->x_max) {
					tbr->x_max = tbr->leafs[j].x;
				}
				if (tbr->leafs[j].y < tbr->y_min)
					tbr->y_min = tbr->leafs[j].y;
				else if (tbr->leafs[j].y > tbr->y_max) {
					tbr->y_max = tbr->leafs[j].y;
				}

				// length
				coord cx = tbr->leafs[j - 1].x - tbr->leafs[j].x;
				coord cy = tbr->leafs[j - 1].y - tbr->leafs[j].y;
				tbr->length[j - 1] = sqrt(cx * cx + cy * cy);
			}

			// boundary
			tbr->offset = (indexer*)malloc(sizeof(indexer) * tbr->count_shapes);
			tbr->xsh_min = (coord*)malloc(sizeof(coord) * tbr->count_shapes);
			tbr->xsh_max = (coord*)malloc(sizeof(coord) * tbr->count_shapes);
			tbr->ysh_min = (coord*)malloc(sizeof(coord) * tbr->count_shapes);
			tbr->ysh_max = (coord*)malloc(sizeof(coord) * tbr->count_shapes);
			indexer k = 0;
			indexer tj = tbr->leafs[0].number;
			tbr->offset[k] = 0;
			coord xsh_min, xsh_max, ysh_min, ysh_max;
			xsh_min = xsh_max = tbr->leafs[0].x;
			ysh_min = ysh_max = tbr->leafs[0].y;
			for (indexer j = 1; j < tbr->count_leafs; ++j) {
				if (tj != tbr->leafs[j].number) {
					tbr->xsh_min[k] = xsh_min;
					tbr->xsh_max[k] = xsh_max;
					tbr->ysh_min[k] = ysh_min;
					tbr->ysh_max[k] = ysh_max;
					++k;
					tbr->offset[k] = j;
					xsh_min = xsh_max = tbr->leafs[j].x;
					ysh_min = ysh_max = tbr->leafs[j].y;
					tj = tbr->leafs[j].number;
					continue;
				}
				if (xsh_min > tbr->leafs[j].x)
					xsh_min = tbr->leafs[j].x;
				else if (xsh_max < tbr->leafs[j].x)
					xsh_max = tbr->leafs[j].x;
				if (ysh_min > tbr->leafs[j].y)
					ysh_min = tbr->leafs[j].y;
				else if (ysh_max < tbr->leafs[j].y)
					ysh_max = tbr->leafs[j].y;
			}
			tbr->xsh_min[k] = xsh_min;
			tbr->xsh_max[k] = xsh_max;
			tbr->ysh_min[k] = ysh_min;
			tbr->ysh_max[k] = ysh_max;
			++k;
			printf("SHAPES %u = %u (%u)\n", i, tbr->count_shapes, k);
#ifdef PRINT_SVG
			//char ch[1024];
			//sprintf_s(ch, 1024, "%u: min X = %f, min Y = %f, max X = %f, max Y = %f", i, tbr->x_max, tbr->y_min, tbr->x_max, tbr->y_max);
			//lprintf(ch);
#endif //PRINT_SVG
		}
		print_file_svg(m_nodes, 1, 0, "7_1.svg");

		//find_branch_centers(m_nodes);
		separate_branches(m_nodes);
		char ch1[64];
#ifndef _WIN
		snprintf(ch1, 64, "7_9.svg");
#else
		sprintf_s(ch1, 64, "7_9.svg");
#endif
		print_file_svg(m_nodes, 1, 0, ch1);
		for (indexer i = 0; i < m_nodes->count_child_nodes; ++i) {
			// char ch1[64];
			//continue;
#ifndef _WIN
			snprintf(ch1, 64, "7_%u.svg", i + 10);
#else
			sprintf_s(ch1, 64, "7_%u.svg", i + 10);
#endif
			print_file_svg((struct node*)(m_nodes->child_node[i]), 1, 0, ch1);
		}

		// finish line
		separate_nodes(m_nodes);
#ifdef PRINT_SVG
/*		for (indexer j = 0; j < m_nodes->count_child_nodes; ++j) {
			struct node *tnd = (struct node*)(m_nodes->child_node[j]);
			sprintf_s(ch1, 64, "7_%u.svg", j + 100);
			print_file_svg(tnd, 1, 0, ch1);

			for (indexer i = 0; i < tnd->count_child_nodes; ++i) {
				// char ch1[64];
				//continue;
				sprintf_s(ch1, 64, "7_%u.svg", j * 100 + i + 1000);
				print_file_svg((struct node*)(tnd->child_node[i]), 1, 0, ch1);
			}
		}*/
#endif // PRINT_SVG

		// testing
		lprintf("start");
		try_find(m_nodes, lll, count_of_leafs);
		lprintf("end");
	}

	free(nds);
	if (lll)
		free(lll);
	free(offsets_leafs);
	del_root();
	lprintf("Done");

#ifdef _WIN
	_getch();
#endif
	return 0;
/*
	// init tree
	init_root(&(lll[0]));

	unsigned ii = 0;
	unsigned offset = 0; // offsets_leafs[ii++];
	// create tree
	for (unsigned i = 0; i < count_of_leafs; i += offset) {
		//printf("=============================== ADD  %u\n", i);
		offset = offsets_leafs[ii++];
		add(m_nodes, true, &(lll[i]), offset);
	}

	free(offsets_leafs);

	//char tch[1024] = {0};
	sprintf_s(tch, 1024, "count = %d, x_min = %f, x_max = %f, y_min = %f, y_max = %f", m_nodes->count_child_nodes, m_nodes->x1, m_nodes->x2, m_nodes->y1, m_nodes->y2);
	lprintf(tch);
*/
/*	FILE *f2;
	errno_t t = fopen_s(&f2, "c:/projects/tmp/1/6_.svg", "w");
	unsigned t1 = 0;

	fprintf(f2, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" height=\"360000px\"  width=\"180000px\">\n"); //  height=\"400px\"  width=\"400px\"
	for (unsigned i = 0; i < count_of_leafs; ++i) {
		fprintf(f2, "\t<polygon points=\"");
		t1 = lll[i].number;
		for (; ; ++i) {
			if (t1 == lll[i].number)
				fprintf(f2, "%u,%u ", (unsigned)lll[i].x, (unsigned)lll[i].y);
			else
				break;
			//printf("%u,%u ", (unsigned)dt[i].pts[j].x, (unsigned)dt[i].pts[j].y);
		}
		//printf("\n");
		fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(0, 0, 0)\"/>\n");
		--i;
	}

	struct node *nd = m_nodes;
	for (unsigned i = 0; i < m_nodes->count_child_nodes; ++i) {
		if (nd->is_last_node) {
			//struct branch *br = (struct branch*)(nd->child_node)[i];
			fprintf(f2, "\t<polygon points=\"");
			fprintf(f2, "%u,%u %u,%u %u,%u %u,%u", (unsigned)(nd->x1), (unsigned)(nd->y1), (unsigned)(nd->x2), (unsigned)(nd->y1), (unsigned)(nd->x2), (unsigned)(nd->y2), (unsigned)(nd->x1), (unsigned)(nd->y2));
			fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(255, 100, 100)\" fill=\"none\"/>\n");
		}
		else {
			// TO DO
		}
	}

	fprintf(f2, "</svg>");
	fclose(f2);
*/
/*	if (lll)
		free(lll);
	//if (m_nodes)
	//	free(m_nodes);
	del_root();

	lprintf("Done");
#ifdef _WIN
	_getch();
#endif
	return 0;
	*/
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

/// init empty root of tree
void init_root(const void* p) {
	m_nodes = (struct node*)malloc(sizeof(struct node));
	m_nodes->count_child_nodes = 1;
	struct branch *tbr = (struct branch*)malloc(sizeof(struct branch));
//#ifndef _WIN
	//*(m_nodes->child_node) = malloc(sizeof(struct branch*) * 1);
//#else
	m_nodes->child_node = (void**)malloc(sizeof(void*) * 1);
//#endif
	tbr->x_min = 0.0;
	tbr->x_max = 0.0;
	tbr->y_min = 0.0;
	tbr->y_max = 0.0;
	(m_nodes->child_node)[0] = tbr;
	/*m_nodes->is_node = (bool*)malloc(sizeof(bool) * 1);
	m_nodes->is_node[0] = false;
	*/
	m_nodes->is_last_node = true;
	m_nodes->center_child_node = (struct center_node_st*)malloc(sizeof(struct center_node_st) * 1);
	m_nodes->center_child_node[0].cy = 0.0;
	m_nodes->center_child_node[0].cx = 0.0;
	m_nodes->center_child_node[0].pos = NULL;

	// set first leaf
	struct leaf *lf = (struct leaf*)p;

	m_nodes->x1 = lf->x;
	m_nodes->y1 = lf->y;
	m_nodes->x2 = lf->x;
	m_nodes->y2 = lf->y;

	//struct branch *br = (struct branch*)m_nodes->child_node[0];
	tbr->count_leafs = 0;
	//br->count_merged_leafs = NULL;
	tbr->merge_next_leaf = NULL;
	tbr->curr_mem_pos = 0;
	tbr->alloc_mem_times = 0;
	tbr->leafs =  NULL;
}

/// init empty root of tree v2
void init_root2(struct node *nd, const void* p) {
	// nd = (struct node*)malloc(sizeof(struct node));
	nd->count_child_nodes = 1;
	struct branch *tbr = (struct branch*)malloc(sizeof(struct branch));
//#ifndef _WIN
	//*(nd->child_node) = malloc(sizeof(struct branch*) * 1);
//#else
	nd->child_node = (void**)malloc(sizeof(void*) * 1);
//#endif
	tbr->x_min = 0.0;
	tbr->x_max = 0.0;
	tbr->y_min = 0.0;
	tbr->y_max = 0.0;
	(nd->child_node)[0] = tbr;
	/*m_nodes->is_node = (bool*)malloc(sizeof(bool) * 1);
	m_nodes->is_node[0] = false;
	*/
	nd->is_last_node = true;
	nd->center_child_node = (struct center_node_st*)malloc(sizeof(struct center_node_st) * 1);
	nd->center_child_node[0].cy = 0.0;
	nd->center_child_node[0].cx = 0.0;
	nd->center_child_node[0].pos = NULL;

	// set first leaf
	struct leaf *lf = (struct leaf*)p;

	nd->x1 = lf->x;
	nd->y1 = lf->y;
	nd->x2 = lf->x;
	nd->y2 = lf->y;

	//struct branch *br = (struct branch*)m_nodes->child_node[0];
	tbr->count_leafs = 0;
	//br->count_merged_leafs = NULL;
	tbr->merge_next_leaf = NULL;
	tbr->curr_mem_pos = 0;
	tbr->alloc_mem_times = 0;
	tbr->leafs = NULL;
}

/// delete root
void del_root()
{
	//const size_t mem_size = 128;
	//unsigned mem_offset = 1;
	//struct node **stack_node = (struct node**)malloc(sizeof(struct node*) * mem_size * mem_offset);
	//unsigned stack_pos = 0;
	//indexer *stack_idx = (indexer*)malloc(sizeof(indexer) * mem_size * mem_offset);

	struct node *nd = m_nodes;

	if (nd->is_last_node) {
		for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
			struct branch *br = (struct branch*)(nd->child_node)[i];
			if (br->merge_next_leaf)
				free(br->merge_next_leaf);
			if (br->leafs)
				free(br->leafs);
			if (br->length)
				free(br->length);
		}
		free((struct branch*)nd->child_node[0]);
		free(nd->child_node);
		free(nd);
		//free(m_nodes->child_node);
		//free(m_nodes);
	}
	else {
		//struct node* nd1 = NULL;
		indexer i = 0;
		struct branch *first_branch = NULL;
		struct node *stack_node[64];
		int stack_pos = 0;
		indexer stack_idx[64];

		while (i < nd->count_child_nodes) {
			if (nd->is_last_node) {
				for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
					struct branch *br = (struct branch*)(nd->child_node)[j];
					if (!first_branch)
						first_branch = br;
					free(br->merge_next_leaf);
					//free(br->center);
					free(br->leafs);
					free(br->length);
					// boundary
					free(br->xsh_max);
					free(br->xsh_min);
					free(br->ysh_max);
					free(br->ysh_min);
				}
				free(nd->child_node);
				free(nd->center_child_node);
				//free(nd);
				// return from stack
				while (stack_pos > 0) {
					stack_pos--;
					nd = stack_node[stack_pos];
					i = stack_idx[stack_pos] + 1;

					if (i < nd->count_child_nodes) {
						stack_idx[stack_pos] = i;
						stack_node[stack_pos] = nd;
						stack_pos++;
						nd = (struct node*)nd->child_node[i];
						i = 0;
						break;
					} else {
						//free(nd->child_node[0]);
						//free(nd->child_node);
						//free(nd->center_child_node);
						//free(nd);
					}
				}
			} else {
				stack_idx[stack_pos] = i;
				stack_node[stack_pos] = nd;
				stack_pos++;
				free(nd->center_child_node);
				nd = (struct node*)nd->child_node[i];
			}
		}
//		return;
		free(first_branch);

		//i = 0;
		nd = m_nodes;
/*		i = 0;
		while (i < nd->count_child_nodes) {
			if (!nd->is_last_node) {
				stack_node[stack_pos] = nd;
				stack_idx[stack_pos] = i;
				stack_pos++;
				i = 0;
				nd = (struct node*)nd->child_node[i];
			} else {
				free(nd->child_node);
				// return from stack
				while (stack_pos > 0) {
					stack_pos--;
					nd = stack_node[stack_pos];
					i = stack_idx[stack_pos] + 1;

					if (i < nd->count_child_nodes) {
						stack_idx[stack_pos] = i;
						stack_node[stack_pos] = nd;
						stack_pos++;
						nd = (struct node*)nd->child_node[i];
						i = 0;
						break;
					}
					free(nd->child_node);
				}
			}
		}
		return;*/

		do /*while (!nd->is_last_node)*/ {
			//stack_idx[stack_pos] = i;
			stack_node[stack_pos] = nd;
			stack_pos++;
			if (!nd->is_last_node)
				nd = (struct node*)nd->child_node[0];
			else {
				//stack_node[stack_pos] = (struct node*)nd->child_node[0];
				//stack_pos++;
				break;
			}
		} while (true);
		for (int i = stack_pos - 1; i >= 0; --i) {
			//free(stack_node[i]->center_child_node);
			free(stack_node[i]);
		}
	}
		///////////////////////////////
	/*	bool flag_del_branches = false;
		do {
			if (i >= nd->count_child_nodes) {
				if (stack_pos > 0) {
					// return from stack
					stack_pos--;
					nd = stack_node[stack_pos];
					i = stack_idx[stack_pos] + 1;
					continue;
				} else {
					break;
				}
			}
			stack_idx[stack_pos] = i;
			stack_node[stack_pos] = nd;
			stack_pos++;
			nd = (struct node*)nd->child_node[i];
			if (nd->is_last_node) {
				struct branch *br;
				for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
					br = (struct branch*)(nd->child_node)[j];
					if (br->merge_next_leaf)
						free(br->merge_next_leaf);
					if (br->leafs)
						free(br->leafs);
				}
				flag_del_branches = true;
				if (!i)
					first_branch = (struct branch*)(nd->child_node)[0];

				free(nd->child_node);
				// return from stack
				stack_pos--;
				nd = stack_node[stack_pos];
				i = stack_idx[stack_pos] + 1;
			} 
		} while (true);
	///	if (first_branch)
	///		free(first_branch);
		free(nd);
	}
	*/
	
	//free(stack_idx);
	//free(stack_node);
	// free(m_nodes->child_node);
	// free(m_nodes);
}

/// add leaf
void add_leaf(struct branch *br, struct leaf *lf, bool merged_next_leaf)
{
	// temp for count
	unsigned tmp = br->count_leafs;
	// struct branch *br = m_nodes[i].child_node;
	if (br->leafs == NULL) {
		br->leafs = (struct leaf*)malloc(sizeof(struct leaf) * 1);
		//br->count_merged_leafs = (unsigned*)malloc(sizeof(unsigned) * 1);
		br->merge_next_leaf = (bool*)malloc(sizeof(bool) * 1);
		br->count_leafs = 1;
	} else {
		++(br->count_leafs);
		//printf("Now leafs = %u\n", br->count_leafs);
		br->leafs = (struct leaf*)realloc(br->leafs, sizeof(struct leaf) * br->count_leafs);
		//br->count_merged_leafs = (unsigned*)realloc(br->count_merged_leafs, sizeof(unsigned) * br->count);
		br->merge_next_leaf = (bool*)realloc(br->merge_next_leaf, sizeof(bool) * br->count_leafs);
	}
	memcpy(br->leafs + tmp, lf, sizeof(struct leaf));
	(br->merge_next_leaf)[tmp] = merged_next_leaf; // count_of_leafs
}

/// add leafs
void add_leafs(struct branch *br, struct leaf *lf, unsigned count_of_leafs)
{
	for (unsigned j = 0; j < count_of_leafs; ++j) {
		if (j == count_of_leafs - 1)
			add_leaf(br, &(lf[j]), false);
		else
			add_leaf(br, &(lf[j]), true);
	}
}

/// add leafs v2
bool add_leafs2(struct branch *br, struct leaf *lf, unsigned count_of_leafs)
{
	// temp for count
	// unsigned tmp = 0; // br->count_leafs;
	// struct branch *br = m_nodes[i].child_node;
	if (br->leafs == NULL) {
		br->leafs = (struct leaf*)malloc(sizeof(struct leaf) * MAX_ADDED_LEAFS_IN_ITER);
		br->merge_next_leaf = (bool*)malloc(sizeof(bool) * MAX_ADDED_LEAFS_IN_ITER);
		br->alloc_mem_times = 1;
	} else {
		if (br->count_leafs + count_of_leafs > br->alloc_mem_times * MAX_ADDED_LEAFS_IN_ITER) {
			++(br->alloc_mem_times);
			br->leafs = (struct leaf*)realloc(br->leafs, sizeof(struct leaf) * br->alloc_mem_times * MAX_ADDED_LEAFS_IN_ITER);
			br->merge_next_leaf = (bool*)realloc(br->merge_next_leaf, sizeof(bool) * br->alloc_mem_times * MAX_ADDED_LEAFS_IN_ITER);
		}
	}
	//br->count_merged_leafs = (unsigned*)malloc(sizeof(unsigned) * 1);
	//br->merge_next_leaf = (bool*)malloc(sizeof(bool) * count_of_leafs);
	br->count_leafs += count_of_leafs;
	memcpy(br->leafs + br->curr_mem_pos, lf, sizeof(struct leaf) * count_of_leafs);
	unsigned i = br->curr_mem_pos;
	for (; i < br->count_leafs - 1; ++i) {
		(br->merge_next_leaf)[i] = true; // count_of_leafs
	}
	(br->merge_next_leaf)[i] = false;
	br->curr_mem_pos = br->count_leafs;
	return true;
}

/// add item to tree (depricate)
/// must call only with root node
/// return flag is added
bool add(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs) {
	if (is_leaf) {
		if (!add_leaf_in_boundary(nd, is_leaf, p, count_of_leafs))
			add_leaf_out_boundary(nd, is_leaf, p, count_of_leafs);

		// check to separate nodes
		// TO DO
//////////////////////////////////////////////////		check_separate(nd);
		//if (((struct branch *)m_nodes[idx_node].child_node)->count_leafs > MAX_ITEMS_IN_NODE) {
			// separate
		//}

		//if (flag)
		//	return true;
		return false;
	} else {
		// not leaf

	}
	return false;
}

/// possible depricated
bool add_leaf_in_boundary(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs)
{
	// flag to added leaf
	//bool flag = false;
	// is leaf
	struct leaf *l = (struct leaf*)p;
	// index of node where insert in boundary
	unsigned idx_node;
	// flag in the bondary
	bool in_boundary = true;

	// recusive to childen nodes
	if (!nd->is_last_node) {
		for  (unsigned i = 0; i < nd->count_child_nodes; ++i) {
			if (add_leaf_in_boundary(&((struct node*)(nd->child_node))[i], is_leaf, p, count_of_leafs))
				return true;
		}
		//lprintf("Error add_leaf_in_boundary 1");
		return false;
	}

	for (unsigned i = 0; i < nd->count_child_nodes; ++i) { /// TO DO may be remove this circle
		/*char tch [1024] = {};
		sprintf(tch, "%d: %f", i, l->x);
		lprintf(tch);
		*/

		// check boundary
		for (unsigned j = 0; j < count_of_leafs; ++j) {
			if (l[j].x >= nd[i].x1 && l[j].x <= nd[i].x2 && l[j].y >= nd[i].y1 && l[j].y <= nd[i].y2) {
				// ???
				idx_node = i;
			} else {
				in_boundary = false;
				break;
			}
		}
		if (in_boundary)
			break;
	}
	
	// debug info
	/* char tch[128] = {0};
	sprintf(tch, "in boundary = %d (x1 = %f, x2 = %f, y1 = %f, y2 = %f)", in_boundary, nd[idx_node].x1, nd[idx_node].x2, nd[idx_node].y1, nd[idx_node].y2);
	lprintf(tch);
	*/

	if (in_boundary) {
		// in the boundary, try to add leaf
		if (!add_leafs2((struct branch *)((nd->child_node)[idx_node]), l, count_of_leafs))
			printf("error 1\n");
		//for (unsigned j = 0; j < count_of_leafs; ++j) {
		//	if (j == count_of_leafs - 1)
		//		add_leaf(m_nodes[i].child_node, &(l[j]), false);
		//	else
		//		add_leaf(m_nodes[i].child_node, &(l[j]), true);
			//struct branch *br = m_nodes[i].child_node;
			/* if (br->leafs == NULL) {
				br->leafs = (struct leaf*)malloc(sizeof(struct leaf) * 1);
				//br->count_merged_leafs = (unsigned*)malloc(sizeof(unsigned) * 1);
				br->merge_next_leaf = (bool*)malloc(sizeof(malloc) * 1);
				br->count = 1;
			} else {
				++(br->count);
				br->leafs = (struct leaf*)realloc(br->leafs, sizeof(struct leaf) * br->count);
				//br->count_merged_leafs = (unsigned*)realloc(br->count_merged_leafs, sizeof(unsigned) * br->count);
				br->merge_next_leaf = (bool*)realloc(br->merge_next_leaf, sizeof(bool) * br->count);
			}
			memcpy(br->leafs + br->count - 1, l, sizeof(struct leaf));
			//br->count_merged_leafs[br->count - 1] = 1; // count_of_leafs
			*/
			//if (count_of_leafs == 1)
			//	br->merge_next_leaf = false;
			//else {
				// TO DO
				// insert next leaf
			//}
		//}
		//flag = true;
		return true;
	}
	//////////// lprintf("Error add_leaf_in_boundary 2");
	return false;
}

/// possible depricated
bool add_leaf_out_boundary(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs)
{
	// TO DO +++
	// not in the boundary
	// unsigned idx_node = 0;
	bool flag = false;
	// is leaf
	struct leaf *l = (struct leaf*)p;

	// need to expand boundary
	coord min_square = (coord)999999999999.0;
	coord min_fake_x1 = (coord)999999999999.0, min_fake_y1 = (coord)999999999999.0, min_fake_x2 = (coord)999999999999.0, min_fake_y2 = (coord)999999999999.0;
	// select node to expand
	//unsigned idx = 0;
	coord fake_x1, fake_y1, fake_x2, fake_y2;

	// as in first time nd == m_nodes than nd->is_last -> true, than nd->child is branch
	//for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
	//	struct node* tnd = (struct node*)((nd->child_node)[j]);
		struct branch* tbr = (struct branch*)((nd->child_node)[0]);
		for (unsigned k = 0; k < count_of_leafs; ++k) {
			if (l[k].x < nd->x1 && l[k].x < nd->x2) {
				// expand left boundary
				fake_x1 = l[k].x;
				if (!k)
					fake_x2 = nd->x2;
			} else if (l[k].x > nd->x1 && l[k].x > nd->x2) {
				// expand right boundary
				if (!k)
					fake_x1 = nd->x1;
				fake_x2 = l[k].x;
			} else if (!k) {
				fake_x1 = nd->x1;
				fake_x2 = nd->x2;
			}
			if (l[k].y < nd->y1 && l[k].y < nd->y2) {
				// expand left boundary
				fake_y1 = l[k].y;
				if (!k)
					fake_y2 = nd->y2;
			} else if (l[k].y > nd->y1 && l[k].y > nd->y2) {
				// expand right boundary
				if (!k)
					fake_y1 = nd->y1;
				fake_y2 = l[k].y;
			} else if (!k) {
				fake_y1 = nd->y1;
				fake_y2 = nd->y2;
			}
		}
		coord square = fabs(fake_x1 - fake_x2) * fabs(fake_y1 - fake_y2);
		char tch[1024] = {0};
///////////		sprintf(tch, "square = %f, x1 = %f, x2 = %f, y1 = %f, y2 = %f, count_of_leafs = %d", square, fake_x1, fake_x2, fake_y1, fake_y2, count_of_leafs);
///////////		lprintf(tch);
		if (square < min_square) {
			// idx_node = j;
			min_square = square;
			min_fake_x1 = fake_x1;
			min_fake_y1 = fake_y1;
			min_fake_x2 = fake_x2;
			min_fake_y2 = fake_y2;
		}
	//}

	// recusive to childen nodes
		// in the first time never work
	if (!nd->is_last_node) {
		//for  (unsigned i = 0; i < nd->count_nodes; ++i) {
			if (add_leaf_out_boundary((struct node*)((nd->child_node)[0/* idx_node */]), is_leaf, p, count_of_leafs)) {
				flag = true;
				//break;
			}
		//}
		//return false;
	} else {
		flag = true;
	}
	
	// TO DO boundary
	if (flag) {
		// expand boundary and add leaf
		// expand node
		// as first time, than never use idx
		nd/*m_nodes[idx_node].*/->x1 = min_fake_x1;
		nd/*m_nodes[idx_node].*/->y1 = min_fake_y1;
		nd/*m_nodes[idx_node].*/->x2 = min_fake_x2;
		nd/*m_nodes[idx_node].*/->y2 = min_fake_y2;
		// add leaf
		if (nd->is_last_node)
			if (!add_leafs2((struct branch *)((nd->child_node)[0/*idx_node*/]), l, count_of_leafs)) {
				printf("error 2\n");
			}
		return true;
	}
	return false;
}

/// old functions (depricated)
struct node* check_separate(struct node* nd)
{
	struct node* tmp_nd = NULL;

	// check leafs
	if (nd->is_last_node) {
		struct branch *br = (struct branch*)nd->child_node;
		//printf("1 = %u > %u\n", br->count_leafs, MAX_ITEMS_IN_NODE);
		if (br != NULL && br->count_leafs > MAX_ITEMS_IN_NODE) {
			// find center
			coord *xc, *yc;
			xc = (coord*)malloc(sizeof(coord) * br->count_leafs);
			yc = (coord*)malloc(sizeof(coord) * br->count_leafs);
			unsigned idxc = 0;
			struct leaf *l = NULL;
			coord minx, maxx, miny, maxy;
			bool first_iter = true;
			//unsigned ixd_center_item = 0;
			//bool find_center_item = false;
			for (unsigned i = 0; i < br->count_leafs; ++i) {
				l = &(br->leafs)[i];
				if (first_iter) {
					minx = l->x;
					miny = l->y;
					maxx = l->x;
					maxy = l->y;
					first_iter = false;
					if ((br->merge_next_leaf)[i])
						continue;
				}

				/// after try to change without if (from habr cuda https://habrahabr.ru/post/204682/)
				if (l->x < minx)
					minx = l->x;
				if (l->x > maxx)
					maxx = l->x;
				if (l->y < miny)
					miny = l->y;
				if (l->y > maxy)
					maxy = l->y;

				if (!(br->merge_next_leaf)[i]) {
					xc[idxc] = (maxx - minx) / 2.0 + minx;
					yc[idxc] = (maxy - miny) / 2.0 + miny;

					// find center in points
					/*if (!find_center_item && br->count_leafs / 2 > i) {
						find_center_item = true;
						ixd_center_item = idxc;
					}*/

					++idxc;
					first_iter = true;
				}
			}
			// try to separate
			// sort xc? than yc
			unsigned count_centers = idxc;
			unsigned *pos_idx_x = (unsigned*)malloc(sizeof(unsigned) * count_centers);
			unsigned *pos_idx_y = (unsigned*)malloc(sizeof(unsigned) * count_centers);
			for (unsigned i = 0; i < count_centers; ++i) {
				pos_idx_x[i] = i;
				pos_idx_y[i] = i;
			}
			// sort xc and yc
			/*quick_sort_r_mass(xc, count_centers, pos_idx_x);
			FILE *f1 = fopen("1.txt", "w");
			if (f1) {
				for (unsigned i = 0; i < count_centers; ++i)
					fprintf(f1, "%d: xc = %f, yc = %f\n", pos_idx_x[i], xc[i], yc[i]);
				fclose(f1);
			} */
///////////////////			q_sort_mas(xc, count_centers, pos_idx_x);
			/* f1 = fopen("2.txt", "w");
			if (f1) {
				for (unsigned i = 0; i < count_centers; ++i)
					fprintf(f1, "%d: xc = %f, yc = %f\n", pos_idx_x[i], xc[i], yc[i]);
				fclose(f1);
			} */
////////////////////			q_sort_mas(yc, count_centers, pos_idx_y);
			// find area axix x
			coord area_x1 = 0.0, area_x2 = 0.0;
			coord area_y1 = 0.0, area_y2 = 0.0;
			area_x1 = (xc[(unsigned)(count_centers / 2)] - xc[0]) * (yc[count_centers - 1] - yc[0]);
			area_x2 = (xc[count_centers - 1] - xc[(unsigned)(count_centers / 2) + 1]) * (yc[count_centers - 1] - yc[0]);
			area_y1 = (yc[(unsigned)(count_centers / 2)] - yc[0]) * (xc[count_centers - 1] - xc[0]);
			area_y2 = (yc[count_centers - 1] - yc[(unsigned)(count_centers / 2) + 1]) * (xc[count_centers - 1] - xc[0]);

			char tch[1024] = {0};
//////////////////////			sprintf(tch, "area_x1 = %f, area_x2 = %f, area_y1 = %f, area_y2 = %f, count_centers = %u", area_x1, area_x2, area_y1, area_y2, count_centers);
			lprintf(tch);
/////////////////////			sprintf(tch, "1: min_x = %f, center1 = %f, center2 = %f, max_x = %f, min_y = %f, max_y = %f", xc[0], xc[(unsigned)(count_centers / 2)], xc[(unsigned)(count_centers / 2) + 1], xc[count_centers - 1], yc[0], yc[count_centers - 1]);
			lprintf(tch);

			// TO DO
			// separate leafs
			// choose separate method (side)
			bool separate_by_x = false;
			coord area_x = area_x1 + area_x2;
			coord area_y = area_y1 + area_y2;
#ifndef _WIN
			min_max_sort(area_x, area_y);
#else
//////////////			min_max_sort(&area_x, &area_y);
#endif
			if (area_x * 1.15 < area_y) {
				if (area_x1 + area_x2 > area_y1 + area_y2) {
					// separate by Y
					separate_by_x = false;
				} else {
					// separate by X
					separate_by_x = true;
				}
			} else {
				// find more square figure
				coord x1 = xc[(unsigned)(count_centers / 2)] - xc[0];
				coord y1 = yc[count_centers - 1] - yc[0];
				coord x2 = xc[count_centers - 1] - xc[0];
				coord y2 = yc[(unsigned)(count_centers / 2)] - yc[0];
#ifndef _WIN
				min_max_sort(x1, y1);
				min_max_sort(x2, y2);
#else
///////////////				min_max_sort(&x1, &y1);
///////////////				min_max_sort(&x2, &y2);
#endif
				if (x1 / y1 > x2 / y2) {
					// separate by Y
					separate_by_x = false;
				} else {
					// separate by X
					separate_by_x = true;
				}
			}
			struct node *tnd = NULL;
			if (separate_by_x)
				tnd = separate_leafs(nd, count_centers, pos_idx_x);
			else
				tnd = separate_leafs(nd, count_centers, pos_idx_y);

			// free
			free(pos_idx_x);
			free(pos_idx_y);
			free(xc);
			free(yc);

			return tnd;
		}
		return NULL;
	} else {
		// recursive moving on tree
		for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
			if (nd->is_last_node)
				break;
			struct node *tnd = (struct node*)(nd->child_node[i]);
			tmp_nd = check_separate(tnd);

			// check temporary node
			if (tmp_nd) {
				// add node to tree
				*(nd->child_node) = realloc(*(nd->child_node), sizeof(struct branch*) * nd->count_child_nodes + 1);
				(nd->child_node)[nd->count_child_nodes] = tmp_nd;
				++(nd->count_child_nodes);
			}
		}

		// check count of nodes
		if (!(nd->is_last_node) && nd->count_child_nodes > MAX_NODES) {
			// separate nodes
			tmp_nd = separate(nd);
			//////////////////////////// to do;
		}
		return tmp_nd;
	}

	return NULL;
}

/// old fubctions (depricated)
struct node* separate_leafs(struct node* nd, unsigned count_senters_items, unsigned *pos_idx) //_x, unsigned *pos_idx_y, bool separate_by_x)
{
	struct branch *br = (struct branch*)nd->child_node;
	struct branch *new_br1 = (struct branch*)malloc(sizeof(struct branch));
	struct branch *new_br2 = (struct branch*)malloc(sizeof(struct branch));

	// check for axis to separate
	/*if (!separate_by_x) {
		unsigned *t = pos_idx_x;
		pos_idx_x = pos_idx_y;
		pos_idx_y = t;
	}*/

	unsigned j;
	for (unsigned i = 0; i < (unsigned)(count_senters_items / 2); ++i) {
		j = 0;
		bool flag = (br->merge_next_leaf)[pos_idx[i] + j];
		do {
			add_leaf(new_br1, br->leafs + pos_idx[i] + j, flag);
			++j;
		} while (flag);
	}
	for (unsigned i = (unsigned)(count_senters_items / 2) + 1; i < count_senters_items; ++i) {
		j = 0;
		bool flag = (br->merge_next_leaf)[pos_idx[i] + j];
		do {
			add_leaf(new_br2, br->leafs + pos_idx[i] + j, flag);
			++j;
		} while (flag);
	}

	// free current branch
	//for (unsigned i = 0; i < br->count_leafs; ++i) {
		free(br->leafs);
		free(br->merge_next_leaf);
	//}
	// update current node
	(nd->child_node)[0] = new_br1;
	// create new node for other branch
	struct node *nd_new = (struct node*)malloc(sizeof(struct node));
	*nd_new->child_node = malloc(sizeof(struct branch*));
	(nd_new->child_node)[0] = new_br2;

	// return new node
	return nd_new;
}

/// compare centers x in branch
int x_cmp(const void* a, const void* b)
{
	struct center_st2 *cena = (struct center_st2*)a;
	struct center_st2 *cenb = (struct center_st2*)b;
	return (int)(cena->cx * 1000 - cenb->cx * 1000);
	return 0;
}

/// compare centers y in branch
int y_cmp(const void* a, const void* b)
{
	struct center_st2 *cena = (struct center_st2*)a;
	struct center_st2 *cenb = (struct center_st2*)b;
	return (int)(cena->cy - cenb->cy);
}

/// separate nodes and leafs into branches
struct node* separate(struct node* nd)
{
	if (!nd->is_last_node)
		return NULL;
	struct branch *br = (struct branch*)(nd->child_node[0]);
	if (br->count_leafs > MAX_ITEMS_IN_NODE) {
		size_t size_separate = (size_t)ceil((double)(br->count_leafs) / MAX_ITEMS_IN_NODE);
		unsigned char t1 = 0, t2 = 0;
		while (size_separate > 0) {
			t1 += size_separate % 2;
			t2++;
			size_separate = size_separate >> 1;
		}
		if (t1 == 1)
			size_separate = (size_t)pow(2, t2 - 1);
		else
			size_separate = (size_t)pow(2, t2);
		printf("Max size = %zu\n", size_separate);

		indexer cleafs = br->count_leafs;
		indexer cshapes = br->count_shapes;
		struct center_st2 *center_start = br->center;
		indexer pos_start = 0;
		struct boudary *positions1 = (struct boudary*)malloc(sizeof(struct boudary) * size_separate);
		struct boudary *positions2 = (struct boudary*)malloc(sizeof(struct boudary) * size_separate);
		unsigned tcleafs = cleafs;

		// correct for easy circle
		//size_separate--;

		positions2[0].idx1 = pos_start;
		positions2[0].idx2 = cshapes - 1;
		indexer j;
		indexer k;
		indexer ccount;

		bool pos_pos = false; // false - used positions1, true = used positions2

		for (unsigned i = 2; i <= size_separate; i = i << 2) {
			tcleafs = cleafs / i;
			pos_pos = false;
			for (k = 0; k < i >> 1; ++k) {
				pos_start = positions2[k].idx1;
				center_start = &(br->center[pos_start]);
				cshapes = positions2[k].idx2 - positions2[k].idx1 + 1;
				ccount = 0;

				// sort shapes by x
				//qsort(br->center, br->count_shapes, sizeof(struct center_st2), x_cmp);
				//qsort_centers_x(br->center, br->count_shapes);
				qsort_centers_x(center_start, cshapes);
				//size_separate--;

				//unsigned tcleafs = cleafs / i;
				//indexer j;
				//indexer k;
				//indexer ccount = 0;
				//for (k = 0; k < i - 1; ++k) {
					ccount = 0;
					for (j = 0; j < cshapes && ccount < tcleafs; ++j) {
						ccount += center_start[j].count_leafs;
					}
					positions1[k * 2].idx1 = pos_start;
					positions1[k * 2].idx2 = pos_start + j - 1;
					positions1[k * 2].count_leafs = ccount;
				// ?	pos_start = j + 1;
				//}
				positions1[k * 2 + 1].idx1 = pos_start + j;
				positions1[k * 2 + 1].idx2 = /* pos_start + j + 1 + */ /* pos_start + cshapes; / == */ positions2[k].idx2;
				ccount = 0;
				for (; j < cshapes; ++j) {
					ccount += center_start[j].count_leafs;
				}
				positions1[k * 2 + 1].count_leafs = ccount;
			}


			// sort shapes by y
			//i = i << 1;
			if (i * 2 > size_separate) {
				break;
			}
			pos_pos = true;
			tcleafs = cleafs / (i * 2);
			for (k = 0; k < i; ++k) {
				pos_start = positions1[k].idx1;
				center_start = &(br->center[pos_start]);
				cshapes = positions1[k].idx2 - positions1[k].idx1 + 1;
				ccount = 0;
				qsort_centers_y(center_start, cshapes);
				// indexer ccount = 0;
				for (j = 0; j < cshapes && ccount < tcleafs; ++j) {
					ccount += center_start[j].count_leafs;
				}
				positions2[k * 2].idx1 = pos_start;
				positions2[k * 2].idx2 = j + pos_start - 1;
				positions2[k * 2].count_leafs = ccount;
				positions2[k * 2 + 1].idx1 = j + pos_start;
				positions2[k * 2 + 1].idx2 = positions1[k].idx2;
				ccount = 0;
				for (; j < cshapes; ++j) {
					ccount += center_start[j].count_leafs;
				}
				positions2[k * 2 + 1].count_leafs = ccount;
			}
			//positions2[k].idx1 = j + 1;

			/*FILE *f1 = NULL;
			fopen_s(&f1, "c:/projects/tmp/1/sort1.txt", "w");
			for (unsigned i = 0; i < br->count_shapes; ++i) {
				fprintf(f1, "%u: x = %f, y = %f\n", i, br->center[i].cx, 0.0); // br->center[i].cy);
			}
			fclose(f1);
			*/
		}

		// copy positions
		if (pos_pos)
			memcpy(positions1, positions2, sizeof(struct boudary) * size_separate);

		// move shapes to other (sorted) positions
		struct branch *br1 = (struct branch*)malloc(sizeof(struct branch) * size_separate);
		indexer total_count_shapes = 0;
		for (unsigned i = 0; i < size_separate; ++i) {
			pos_start = positions1[i].idx1;
			//center_start = &(center_start[pos_start]);
			cshapes = positions1[i].idx2 - positions1[i].idx1 + 1;
			tcleafs = positions1[i].count_leafs;

			br1[i].leafs = (struct leaf*)malloc(sizeof(struct leaf) * tcleafs);
			br1[i].merge_next_leaf = (bool*)malloc(sizeof(bool) * tcleafs);

			br1[i].alloc_mem_times = (unsigned)-1;
			br1[i].curr_mem_pos = (unsigned)-1;
			br1[i].count_shapes = cshapes;
			br1[i].count_leafs = tcleafs;
			indexer count_of_leafs = 0;
			for (indexer j = 0; j < cshapes; ++j) {
				memcpy(br1[i].leafs + count_of_leafs, br->center[j + total_count_shapes].pos_leaf, sizeof(struct leaf) * br->center[j + total_count_shapes].count_leafs);
				indexer k = count_of_leafs;
				for (; k < count_of_leafs + br->center[j + total_count_shapes].count_leafs - 1; ++k) {
					br1[i].merge_next_leaf[k] = true;
				}
				if (k < br1[i].count_leafs) // ???
					br1[i].merge_next_leaf[k] = false;
				count_of_leafs += br->center[j + total_count_shapes].count_leafs;
				// br->center[j + total_count_shapes].moved = true;
			}
			total_count_shapes += cshapes;

			// calculate boundary of branches
			coord xmi = br1[i].leafs[0].x;
			coord xma = xmi;
			coord ymi = br1[i].leafs[0].y;
			coord yma = ymi;
			for (indexer j = 1; j < br1[i].count_leafs; ++j) {
				if (br1[i].leafs[j].x < xmi)
					xmi = br1[i].leafs[j].x;
				if (br1[i].leafs[j].x > xma)
					xma = br1[i].leafs[j].x;
				if (br1[i].leafs[j].y < ymi)
					ymi = br1[i].leafs[j].y;
				if (br1[i].leafs[j].y > yma)
					yma = br1[i].leafs[j].y;
			}
			br1[i].x_min = xmi;
			br1[i].x_max = xma;
			br1[i].y_min = ymi;
			br1[i].y_max = yma;
		}

		// delete old branch and exchange on curernt branches
		free(br->center);
		free(br->leafs);
		free(br->merge_next_leaf);
		free(br);
		free(nd->child_node);
		nd->child_node = (void**)malloc(sizeof(void*) * size_separate);
		free(nd->center_child_node); // ???
		nd->center_child_node = (struct center_node_st*)malloc(sizeof(struct center_node_st) * size_separate);
		for (indexer i = 0; i < size_separate; ++i) {
			nd->child_node[i] = &(br1[i]);
			// centers
			nd->center_child_node[i].cx = (br1[i].x_max - br1[i].x_min) / (coord)2.0 + br1[i].x_min;
			nd->center_child_node[i].cy = (br1[i].y_max - br1[i].y_min) / (coord)2.0 + br1[i].y_min;
			nd->center_child_node[i].pos = &(br1[i]);
		}
		nd->count_child_nodes = (unsigned)size_separate;

		free(positions1);
		free(positions2);

		return nd;
	}

	return NULL;
}

/// create threads for first node
bool create_first_thread(struct node* nd, struct leaf* leafs, unsigned *offsets_leafs, unsigned count_cpus, unsigned count_leafs)
{
	// http://learnc.info/c/pthreads_create_and_join.html
	// https://msdn.microsoft.com/en-us/library/kdzttdcb.aspx
	// pthread_create for linux, _beginthread for windows

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
	uintptr_t *ptr1 = (uintptr_t*)malloc(sizeof(uintptr_t) * count_cpus);
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
		ptr1[i] = _beginthread(first_thread_v2, 0, &(ftst[i]));
	}
	for (unsigned i = 0; i < count_cpus; ++i) {
		WaitForSingleObject((HANDLE)ptr1[i], INFINITE);
	}

	free(idx);
	free(ftst);
	free(ptr1);
#endif // _WIN
	return false;
}

/// firts thread, possible depricated
void first_thread(void *params)
{
#ifndef _WIN
#else

	/* // test
	printf("%u\n", *(unsigned*)params);
	Sleep(1000);
	*/
#endif // _WIN

	struct first_thr_st *ftst = (struct first_thr_st*)params;
	
	char tch[32];
#ifndef _WIN
	snprintf(tch, 32, "Start %u", ftst->start_pos_leafs);
#else
	sprintf_s(tch, 32, "Start %u", ftst->start_pos_leafs);
#endif
	lprintf(tch);

	// init tree
	init_root2(ftst->node_, &(ftst->leafs_)[ftst->start_pos_leafs]);

	unsigned ii = 0;
	unsigned offset = 0; // offsets_leafs[ii++];
	unsigned t1 = 0;
						 // create tree
	for (unsigned i = ftst->start_pos_leafs /*0*/; i < ftst->start_pos_leafs + ftst->count_leaf /*count_of_leafs*/; i += offset) {
		//printf("=============================== ADD  %u\n", i);
		offset = (ftst->offsets_leafs_)[ii++]; // offsets_leafs[ii++];
		add(ftst->node_, true, &((ftst->leafs_)[i]), offset);
		t1 += offset;
	}

#ifndef _WIN
	snprintf(tch, 32, "Stop %u (t1 = %u)", ftst->start_pos_leafs, t1);
#else
	sprintf_s(tch, 32, "Stop %u (t1 = %u)", ftst->start_pos_leafs, t1);
#endif
	lprintf(tch);
	return;
}

/// first tharead v2
#ifndef _WIN
void* first_thread_v2(void *params)
#else
void first_thread_v2(void *params)
#endif
{
#ifndef _WIN
#else

	/* // test
	printf("%u\n", *(unsigned*)params);
	Sleep(1000);
	*/
#endif // _WIN

	struct first_thr_st *ftst = (struct first_thr_st*)params;

	char tch[32];
#ifndef _WIN
	snprintf(tch, 32, "Start %u", ftst->start_pos_leafs);
#else
	sprintf_s(tch, 32, "Start %u", ftst->start_pos_leafs);
#endif
	lprintf(tch);

	// init tree
	init_root2(ftst->node_, &(ftst->leafs_)[ftst->start_pos_leafs]);

	unsigned ii = 0;
	unsigned offset = 0; // offsets_leafs[ii++];
	unsigned t1 = 0;
	// create tree
	for (unsigned i = ftst->start_pos_leafs /*0*/; i < ftst->start_pos_leafs + ftst->count_leaf /*count_of_leafs*/; i += offset) {
		//printf("=============================== ADD  %u\n", i);
		offset = (ftst->offsets_leafs_)[ii++]; // offsets_leafs[ii++];
		// add(ftst->node_, true, &((ftst->leafs_)[i]), offset);
		// if (!add_leafs2((struct branch *)((nd->child_node)[idx_node]), l, count_of_leafs))
		if (!add_leafs2((struct branch *)((ftst->node_->child_node)[0]), &((ftst->leafs_)[i]), offset))
			printf("error 1\n");
		t1 += offset;
	}

#ifndef _WIN
	snprintf(tch, 32, "Stop %u (t1 = %u)", ftst->start_pos_leafs, t1);
	lprintf(tch);
	pthread_exit(0);
	return NULL;

#else
	sprintf_s(tch, 32, "Stop %u (t1 = %u)", ftst->start_pos_leafs, t1);
	lprintf(tch);
	return;
#endif
}

/// print small data in svg file
void print_file_svg(struct node *nds, unsigned count_nodes, unsigned count_shapes, const char *file_name)
{
#ifdef PRINT_SVG
	FILE *f2;
	char name[256];
	sprintf_s(name, 256, "c:/projects/tmp/1/%s", file_name);
	errno_t t = fopen_s(&f2, name, "w");
	unsigned t1 = 0;

	fprintf(f2, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" height=\"10000px\"  width=\"10000px\">\n"); //  height=\"400px\"  width=\"400px\"
	for (unsigned k = 0; k < count_nodes; ++k) {
		struct node *nd = &(nds[k]);
		if (!nd->is_last_node) {
			// draw boundary
			//break;
			fprintf(f2, "\t<polygon points=\"");
			fprintf(f2, "%u,%u %u,%u %u,%u %u,%u", (unsigned)(nd->x1), (unsigned)(nd->y1), (unsigned)(nd->x2), (unsigned)(nd->y1), (unsigned)(nd->x2), (unsigned)(nd->y2), (unsigned)(nd->x1), (unsigned)(nd->y2));
			fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(200, 55, 75)\" fill=\"none\"/>\n"); // rgb(255, 55, 100)
		}
		else {
			for (indexer j = 0; j < nd->count_child_nodes; ++j) {
				struct branch *br = (struct branch*)(nd->child_node)[j];
				char color[64];
				sprintf_s(color, 64, "rgb(%u,%u,%u)", rand() % 100 + 155, rand() % 100 + 155, rand() % 100 + 100);

				for (unsigned i = 0; i < br->count_leafs /*count_of_leafs*/; ++i) {
					fprintf(f2, "\t<polygon points=\"");
					t1 = (br->leafs)[i].number; // lll[i].number;
					for (; i < br->count_leafs; ++i) {
						if (t1 == (br->leafs)[i].number /*lll[i].number*/)
							fprintf(f2, "%u,%u ", (unsigned)(br->leafs)[i].x /*lll[i].x*/, (unsigned)(br->leafs)[i].y /*lll[i].y*/);
						else
							break;
						//printf("%.02f,%.02f ", /*(unsigned)*/(br->leafs)[i].x, /*(unsigned)*/(br->leafs)[i].y);
					}
					//printf("\n");
					fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(0, 0, 0)\" fill=\"%s\"/>\n", color); // rgb(0, 0, 0) rgb(150,150,255)
					--i;
				}

				// draw centers
				if (count_shapes) {
					for (unsigned i = 0; i < count_shapes; ++i) {
						fprintf(f2, "\t<circle cx=\"%u\" cy=\"%u\" r=\"1\" stroke-width=\"1\" stroke=\"rgb(50, 255, 100)\"/>\n", (unsigned)(br->center[i].cx), (unsigned)(br->center[i].cy));
					}
				}
			}

			// node boundary
			fprintf(f2, "\t<polygon points=\"");
			fprintf(f2, "%u,%u %u,%u %u,%u %u,%u", (unsigned)(nd->x1), (unsigned)(nd->y1), (unsigned)(nd->x2), (unsigned)(nd->y1), (unsigned)(nd->x2), (unsigned)(nd->y2), (unsigned)(nd->x1), (unsigned)(nd->y2));
			fprintf(f2, "\" stroke-width=\"2\" stroke=\"rgb(255, 100, 255)\" fill=\"none\"/>\n");

			// struct node *nd = m_nodes;
			if (nd->is_last_node) {
				for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
					struct branch *br = (struct branch*)(nd->child_node)[i];
					// draw boundary branch
					fprintf(f2, "\t<polygon points=\"");
					fprintf(f2, "%u,%u %u,%u %u,%u %u,%u", (unsigned)(br->x_min), (unsigned)(br->y_min), (unsigned)(br->x_max), (unsigned)(br->y_min), (unsigned)(br->x_max), (unsigned)(br->y_max), (unsigned)(br->x_min), (unsigned)(br->y_max));
					fprintf(f2, "\" stroke-width=\"1\" stroke=\"rgb(200, 55, 75)\" fill=\"none\"/>\n"); // rgb(255, 55, 100)
				}
			}
			else {
				// TO DO
			}
		}
	}

	fprintf(f2, "</svg>");
	fclose(f2);
#endif //PRINT_SVG
}

/// find centers (nodes, leafs in branches)
bool find_centers(const struct node *nd, const unsigned count_shapes)
{
	if (!nd->is_last_node) {
		// current item is node
		return false;
	}

	coord tx, ty;
	struct branch* br = (struct branch*)(nd->child_node[0]);

	br->center = (struct center_st2*)malloc(sizeof(struct center_st2) * count_shapes);
	// hope that this is the fitst time and alloc memory
	//////////////////br->count_merged_pos_leafs = (unsigned*)malloc(sizeof(unsigned) * count_shapes);
	// br->center.pos_leaf = (struct leaf**)malloc(sizeof(struct leaf*) * count_shapes);
	//br->center->cx = (coord*)malloc(sizeof(coord) * count_shapes);
	//br->center->cy = (coord*)malloc(sizeof(coord) * count_shapes);

	unsigned j = 0; // count leafs in union
	unsigned k = 0; // index of shape
	br->center[k].pos_leaf = &(br->leafs[0]);
	//j++;
	for (unsigned i = 0; i < br->count_leafs; ++i) {
		if (!br->merge_next_leaf[i]) {
			// calculate center
			tx = 0.0;
			ty = 0.0;
			for (unsigned t1 = i - j; t1 <= i; ++t1) {
				tx += br->leafs[t1].x;
				ty += br->leafs[t1].y;
			}
			j++;
			//br->count_merged_pos_leafs[k] = j;
			br->center[k].cx = tx / j;
			br->center[k].cy = ty / j;
			br->center[k].count_leafs = j;
			// br->center[k].moved = false;
			
			// go to the next shape
			k++;
			if (k >= count_shapes)
				break;
			j = 0;
			// i++;
			br->center[k].pos_leaf = &(br->leafs[i + 1]);
		} else {
			j++;
		}
	}
	br->count_shapes = count_shapes;
	return true;
}

/// calculate branch centers (depricated)
bool find_branch_centers(struct node *nd)
{
	if (!nd->is_last_node) {
		// not current step
		return false;
	}

/*	for (indexer i = 0; i < nd->count_child_nodes; ++i) {
		struct branch *br = (struct branch*)nd->child_node[i];
		br->branch_center.cx = br->x_max - br->y_min;
		br->branch_center.cy = br->y_max - br->y_min;
	}
	nd->center.cx = nd->x2 - nd->x1;
	nd->center.cy = nd->y2 - nd->y1;
	*/
	return false;
}

/// separate branches
struct node* separate_branches(struct node* nd)
{
	if (!nd->is_last_node) {
		return NULL;
	}

	if (nd->count_child_nodes <= MAX_NODES) {
		return nd;
	}

	size_t size_separate = (size_t)ceil((double)(nd->count_child_nodes) / MAX_NODES);
	unsigned char t1 = 0, t2 = 0;
	while (size_separate > 0) {
		t1 += size_separate % 2;
		t2++;
		size_separate = size_separate >> 1;
	}
	if (t1 == 1)
		size_separate = (size_t)pow(2, t2 - 1);
	else
		size_separate = (size_t)pow(2, t2);
	printf("Max size branches in node = %zu\n", size_separate);

	indexer cbranches = nd->count_child_nodes;
	struct center_node_st *center_start = nd->center_child_node;
	indexer pos_start = 0;
	struct boudary *positions1 = (struct boudary*)malloc(sizeof(struct boudary) * size_separate);
	struct boudary *positions2 = (struct boudary*)malloc(sizeof(struct boudary) * size_separate);
	unsigned tcbranches = cbranches;
	unsigned brinscope = cbranches;

	positions2[0].idx1 = pos_start;
	positions2[0].idx2 = cbranches - 1;
	//indexer j;
	indexer k;
	//indexer ccount;

	bool pos_pos = false; // false - used positions1, true = used positions2

	for (unsigned i = 2; i <= size_separate; i = i << 2) {
		tcbranches = cbranches / i; // ???
		pos_pos = false;
		for (k = 0; k < i >> 1; ++k) {
			pos_start = positions2[k].idx1;
			center_start = &(nd->center_child_node[pos_start]);
			brinscope = positions2[k].idx2 - positions2[k].idx1 + 1;
			//ccount = 0;

			// sort shapes by x
			//qsort(br->center, br->count_shapes, sizeof(struct center_st2), x_cmp);
			//qsort_centers_x(br->center, br->count_shapes);
			qsort_node_centers_x(center_start, brinscope);
			//size_separate--;
			
			positions1[k * 2].idx1 = pos_start;
			positions1[k * 2].idx2 = pos_start + (brinscope >> 1) - 1;
			positions1[k * 2].count_leafs = (brinscope >> 1);
			// ?	pos_start = j + 1;
			//}
			positions1[k * 2 + 1].idx1 = pos_start + (brinscope >> 1);
			positions1[k * 2 + 1].idx2 = /* pos_start + j + 1 + */ /* pos_start + cshapes; / == */ positions2[k].idx2;
			positions1[k * 2 + 1].count_leafs = brinscope - (brinscope >> 1);
		}


		// sort shapes by y
		//i = i << 1;
		if (i * 2 > size_separate) {
			break;
		}
		pos_pos = true;
		tcbranches = cbranches / (i * 2);  // ???
		for (k = 0; k < i; ++k) {
			pos_start = positions1[k].idx1;
			center_start = &(nd->center_child_node[pos_start]);
			brinscope = positions1[k].idx2 - positions1[k].idx1 + 1;
			//ccount = 0;
			qsort_node_centers_y(center_start, brinscope);

			positions2[k * 2].idx1 = pos_start;
			positions2[k * 2].idx2 = pos_start + (brinscope >> 1) - 1;
			positions2[k * 2].count_leafs = (brinscope >> 1);
			positions2[k * 2 + 1].idx1 = pos_start + (brinscope >> 1);
			positions2[k * 2 + 1].idx2 = positions1[k].idx2;
			positions2[k * 2 + 1].count_leafs = brinscope - (brinscope >> 1);
		}
		//positions2[k].idx1 = j + 1;
	}

	// copy positions
	if (pos_pos)
		memcpy(positions1, positions2, sizeof(struct boudary) * size_separate);
	
	// move branches to other (sorted nodes) positions
	struct branch* tbr1 = (struct branch*)malloc(sizeof(struct branch) * nd->count_child_nodes);
	for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
		memcpy(&(tbr1[i]), nd->center_child_node[i].pos, sizeof(struct branch));
	}

	// prepare new nodes
	struct node *nd1 = (struct node*)malloc(sizeof(struct node) * size_separate);
	//indexer total_count_br = 0;
	for (unsigned i = 0; i < size_separate; ++i) {
		init_root2(&(nd1[i]), &(((struct branch*)(nd->child_node[0]))->leafs[0]));
		// assign branches to node
		free((struct branch*)nd1[i].child_node[0]);
		free(nd1[i].child_node);
		nd1[i].child_node = (void**)malloc(sizeof(void*) * (positions1[i].idx2 - positions1[i].idx1 + 1));
		for (unsigned j = positions1[i].idx1; j < positions1[i].idx2 + 1; ++j) {
			nd1[i].child_node[j - positions1[i].idx1] = &(tbr1[j]);
		}
		nd1[i].count_child_nodes = positions1[i].idx2 + 1 - positions1[i].idx1;
		
		// calculate new boundary
		coord xmi = tbr1[positions1[i].idx1].x_min;
		coord xma = tbr1[positions1[i].idx1].x_max;
		coord ymi = tbr1[positions1[i].idx1].y_min;
		coord yma = tbr1[positions1[i].idx1].y_max;
		for (unsigned j = positions1[i].idx1 + 1; j < positions1[i].idx2 + 1; ++j) {
			if (tbr1[j].x_min < xmi) {
				xmi = tbr1[j].x_min;
			}
			if (tbr1[j].x_max > xma) {
				xma = tbr1[j].x_max;
			}
			if (tbr1[j].y_min < ymi) {
				ymi = tbr1[j].y_min;
			}
			if (tbr1[j].y_max > yma) {
				yma = tbr1[j].y_max;
			}
		}
		nd1[i].x1 = xmi;
		nd1[i].x2 = xma;
		nd1[i].y1 = ymi;
		nd1[i].y2 = yma;
	}
	
	// assign to root
	//free(br);
	//free(nd->child_node);
	//nd->child_node = (void**)malloc(sizeof(void*) * size_separate);

	// free branches and nodes
	free(nd->child_node);
	nd->child_node = (void**)malloc(sizeof(void*) * size_separate);
	// free centers
	free(nd->center_child_node);
	nd->center_child_node = (struct center_node_st*)malloc(sizeof(struct center_node_st) * size_separate);
	for (indexer i = 0; i < size_separate; ++i) {
		nd->child_node[i] = &(nd1[i]);
		// centers
		nd->center_child_node[i].cx = (nd1[i].x2 - nd1[i].x1) / (coord)2.0 + nd1[i].x1;
		nd->center_child_node[i].cy = (nd1[i].y2 - nd1[i].y1) / (coord)2.0 + nd1[i].y1;
		nd->center_child_node[i].pos = &(nd1[i]);
	}
	nd->count_child_nodes = (unsigned)size_separate;
	nd->is_last_node = false;

	free(positions1);
	free(positions2);

	return nd;
////////////////////////////////////////////////////////////////

	return NULL;
}

/// separate nodes
struct node* separate_nodes(struct node* nd)
{
	if (nd->is_last_node) {
		return NULL;
	}

	if (nd->count_child_nodes <= MAX_NODES) {
		return nd;
	}

	while (nd->count_child_nodes > MAX_NODES) {
		size_t size_separate = (size_t)ceil((double)(nd->count_child_nodes) / MAX_NODES);
		unsigned char t1 = 0, t2 = 0;
		while (size_separate > 0) {
			t1 += size_separate % 2;
			t2++;
			size_separate = size_separate >> 1;
		}
		if (t1 == 1)
			size_separate = (size_t)pow(2, t2 - 1);
		else
			size_separate = (size_t)pow(2, t2);
		printf("Max size nodes in node = %zu\n", size_separate);

		indexer cnodes = nd->count_child_nodes;
		struct center_node_st *center_start = nd->center_child_node;
		indexer pos_start = 0;
		struct boudary *positions1 = (struct boudary*)malloc(sizeof(struct boudary) * size_separate);
		struct boudary *positions2 = (struct boudary*)malloc(sizeof(struct boudary) * size_separate);
		unsigned tcnodes = cnodes;
		unsigned ndinscope = cnodes;

		positions2[0].idx1 = pos_start;
		positions2[0].idx2 = cnodes - 1;
		//indexer j;
		indexer k;
		indexer ccount;

		bool pos_pos = false; // false - used positions1, true = used positions2

		for (unsigned i = 2; i <= size_separate; i = i << 2) {
			tcnodes = cnodes / i; // ???
			pos_pos = false;
			for (k = 0; k < i >> 1; ++k) {
				pos_start = positions2[k].idx1;
				center_start = &(nd->center_child_node[pos_start]);
				ndinscope = positions2[k].idx2 - positions2[k].idx1 + 1;
				ccount = 0;

				// sort shapes by x
				//qsort(br->center, br->count_shapes, sizeof(struct center_st2), x_cmp);
				//qsort_centers_x(br->center, br->count_shapes);
				qsort_node_centers_x(center_start, ndinscope);
				//size_separate--;

				positions1[k * 2].idx1 = pos_start;
				positions1[k * 2].idx2 = pos_start + (ndinscope >> 1) - 1;
				positions1[k * 2].count_leafs = (ndinscope >> 1);
				// ?	pos_start = j + 1;
				//}
				positions1[k * 2 + 1].idx1 = pos_start + (ndinscope >> 1);
				positions1[k * 2 + 1].idx2 = /* pos_start + j + 1 + */ /* pos_start + cshapes; / == */ positions2[k].idx2;
				positions1[k * 2 + 1].count_leafs = ndinscope - (ndinscope >> 1);
			}


			// sort shapes by y
			//i = i << 1;
			if (i * 2 > size_separate) {
				break;
			}
			pos_pos = true;
			tcnodes = cnodes / (i * 2);  // ???
			for (k = 0; k < i; ++k) {
				pos_start = positions1[k].idx1;
				center_start = &(nd->center_child_node[pos_start]);
				ndinscope = positions1[k].idx2 - positions1[k].idx1 + 1;
				ccount = 0;
				qsort_node_centers_y(center_start, ndinscope);

				positions2[k * 2].idx1 = pos_start;
				positions2[k * 2].idx2 = pos_start + (ndinscope >> 1) - 1;
				positions2[k * 2].count_leafs = (ndinscope >> 1);
				positions2[k * 2 + 1].idx1 = pos_start + (ndinscope >> 1);
				positions2[k * 2 + 1].idx2 = positions1[k].idx2;
				positions2[k * 2 + 1].count_leafs = ndinscope - (ndinscope >> 1);
			}
			//positions2[k].idx1 = j + 1;
		}

		// copy positions
		if (pos_pos)
			memcpy(positions1, positions2, sizeof(struct boudary) * size_separate);

		// move branches to other (sorted nodes) positions
		struct node* tnd1 = (struct node*)malloc(sizeof(struct node) * nd->count_child_nodes);
		for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
			memcpy(&(tnd1[i]), nd->center_child_node[i].pos, sizeof(struct node));
		}

		// prepare new nodes
		struct node *nd1 = (struct node*)malloc(sizeof(struct node) * size_separate);
		indexer total_count_br = 0;
		for (unsigned i = 0; i < size_separate; ++i) {
			struct leaf lf;
			lf.x = 0.0;
			lf.y = 0.0;
			lf.number = 0;
			init_root2(&(nd1[i]), &lf);
			// assign branches to node
			free((struct branch*)nd1[i].child_node[0]);
			free(nd1[i].child_node);
			nd1[i].child_node = (void**)malloc(sizeof(void*) * (positions1[i].idx2 - positions1[i].idx1 + 1));
			for (unsigned j = positions1[i].idx1; j < positions1[i].idx2 + 1; ++j) {
				nd1[i].child_node[j - positions1[i].idx1] = &(tnd1[j]);
			}
			nd1[i].count_child_nodes = positions1[i].idx2 + 1 - positions1[i].idx1;

			// calculate new boundary
			coord xmi = tnd1[positions1[i].idx1].x1;
			coord xma = tnd1[positions1[i].idx1].x2;
			coord ymi = tnd1[positions1[i].idx1].y1;
			coord yma = tnd1[positions1[i].idx1].y2;
			for (unsigned j = positions1[i].idx1 + 1; j < positions1[i].idx2 + 1; ++j) {
				if (tnd1[j].x1 < xmi) {
					xmi = tnd1[j].x1;
				}
				if (tnd1[j].x2 > xma) {
					xma = tnd1[j].x2;
				}
				if (tnd1[j].y1 < ymi) {
					ymi = tnd1[j].y1;
				}
				if (tnd1[j].y2 > yma) {
					yma = tnd1[j].y2;
				}
			}
			nd1[i].x1 = xmi;
			nd1[i].x2 = xma;
			nd1[i].y1 = ymi;
			nd1[i].y2 = yma;
			nd1[i].is_last_node = false;
		}

		// assign to root
		//free(br);
		//free(nd->child_node);
		//nd->child_node = (void**)malloc(sizeof(void*) * size_separate);

		// free branches and nodes
		free(nd->child_node);
		nd->child_node = (void**)malloc(sizeof(void*) * size_separate);
		// free centers
		free(nd->center_child_node);
		nd->center_child_node = (struct center_node_st*)malloc(sizeof(struct center_node_st) * size_separate);
		for (indexer i = 0; i < size_separate; ++i) {
			nd->child_node[i] = &(nd1[i]);
			// centers
			nd->center_child_node[i].cx = (nd1[i].x2 - nd1[i].x1) / (coord)2.0 + nd1[i].x1;
			nd->center_child_node[i].cy = (nd1[i].y2 - nd1[i].y1) / (coord)2.0 + nd1[i].y1;
			nd->center_child_node[i].pos = &(nd1[i]);
		}
		nd->count_child_nodes = (unsigned)size_separate;
		nd->is_last_node = false;

		free(positions1);
		free(positions2);
	}

	return nd;
	////////////////////////////////////////////////////////////////

	return NULL;
}
