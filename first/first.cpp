#include "pch.h"
#include <stdio.h>
//#include <malloc.h>
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
//#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "unimem.h"
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
unsigned m_max_added_leafs_in_iter = 100000;


/// main function
struct node* create_rtree(struct leaf* lll, unsigned count_of_leafs, unsigned *offsets_leafs, unsigned count_shapes)
{
	// threads
	unsigned cpus = 8;
#ifndef _WIN
	cpus = sysconf(_SC_NPROCESSORS_ONLN);
	char tch[64];
	snprintf(tch, 64, "CPUS = %u", cpus);
#else
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	cpus = sysinfo.dwNumberOfProcessors;
	char tch[64];
	sprintf_s(tch, 64, "CPUS = %u", cpus);
#endif
	cpus = 1;
	lprintf(tch);

	// calculate maximal size for adding leafs in iteration for malloc and realloc
	m_max_added_leafs_in_iter = (unsigned)(count_of_leafs / cpus * 0.3);

	alignas(16) struct node *nds = (struct node*)aligned_alloc(16, sizeof(struct node) * cpus);
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
	total_mem = (unsigned)(ceil((double)total_mem / 16.0) * 16);
	//((struct branch*)(nds[0].child_node[0]))->leafs = (struct leaf*)realloc(((struct branch*)(nds[0].child_node[0]))->leafs, sizeof(struct leaf) * total_mem);
	init_root(&(lll[0]));
	struct branch *br1 = (struct branch*)(m_nodes->child_node[0]);
#ifdef OLD_LEAFS
	br1->leafs = (struct leaf*)aligned_alloc(16, sizeof(struct leaf) * total_mem);
#else
	// TO DO LEAFS
	br1->leaf_x = (coord*)aligned_alloc(16, sizeof(coord) * total_mem);
	br1->leaf_y = (coord*)aligned_alloc(16, sizeof(coord) * total_mem);
	br1->leaf_number = (indexer*)aligned_alloc(16, sizeof(indexer) * total_mem);
#endif // OLD_LEAFS
	br1->merge_next_leaf = (bool*)aligned_alloc(16, sizeof(bool) * total_mem);
	for (unsigned i = 0; i < cpus; ++i) {
		// copy leafs
#ifdef OLD_LEAFS
		memcpy(br1->leafs + br1->count_leafs, ((struct branch*)(nds[i].child_node[0]))->leafs, ((struct branch*)(nds[i].child_node[0]))->count_leafs * sizeof(struct leaf));
#else
		// TO DO LEAFS
		memcpy(br1->leaf_x + br1->count_leafs, ((struct branch*)(nds[i].child_node[0]))->leaf_x, ((struct branch*)(nds[i].child_node[0]))->count_leafs * sizeof(coord));
		memcpy(br1->leaf_y + br1->count_leafs, ((struct branch*)(nds[i].child_node[0]))->leaf_y, ((struct branch*)(nds[i].child_node[0]))->count_leafs * sizeof(coord));
		memcpy(br1->leaf_number + br1->count_leafs, ((struct branch*)(nds[i].child_node[0]))->leaf_number, ((struct branch*)(nds[i].child_node[0]))->count_leafs * sizeof(indexer));
#endif // OLD_LEAFS
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
		_aligned_free(((struct branch*)(nds[i].child_node[0]))->merge_next_leaf);
#ifdef OLD_LEAFS
		_aligned_free(((struct branch*)(nds[i].child_node[0]))->leafs);
#else
		// TO DO LEAFS
		_aligned_free(((struct branch*)(nds[i].child_node[0]))->leaf_x);
		_aligned_free(((struct branch*)(nds[i].child_node[0]))->leaf_y);
		_aligned_free(((struct branch*)(nds[i].child_node[0]))->leaf_number);
#endif // OLD_LEAFS
		_aligned_free(((struct branch*)(nds[i].child_node[0])));
		_aligned_free(nds[i].center_child_node);
		_aligned_free(nds[i].child_node);
	}
	_aligned_free(nds);

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
			tbr->length = (coord*)aligned_alloc(16, sizeof(coord) * tbr->count_leafs);

#ifdef OLD_LEAFS
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
#else
			// TO DO LEAFS
			tbr->x_min = tbr->leaf_x[0];
			tbr->x_max = tbr->leaf_x[0];
			tbr->y_min = tbr->leaf_y[0];
			tbr->y_max = tbr->leaf_y[0];
			for (indexer j = 1; j < tbr->count_leafs; ++j) {
				// boundary
				if (tbr->leaf_x[j] < tbr->x_min)
					tbr->x_min = tbr->leaf_x[j];
				else if (tbr->leaf_x[j] > tbr->x_max) {
					tbr->x_max = tbr->leaf_x[j];
				}
				if (tbr->leaf_y[j] < tbr->y_min)
					tbr->y_min = tbr->leaf_y[j];
				else if (tbr->leaf_y[j] > tbr->y_max) {
					tbr->y_max = tbr->leaf_y[j];
				}

				// length
				coord cx = tbr->leaf_x[j - 1] - tbr->leaf_x[j];
				coord cy = tbr->leaf_y[j - 1] - tbr->leaf_y[j];
				tbr->length[j - 1] = sqrt(cx * cx + cy * cy);
			}
#endif // OLD_LEAFS

			// boundary
			tbr->offset = (indexer*)aligned_alloc(16, sizeof(indexer) * tbr->count_shapes);
			tbr->xsh_min = (coord*)aligned_alloc(16, sizeof(coord) * tbr->count_shapes);
			tbr->xsh_max = (coord*)aligned_alloc(16, sizeof(coord) * tbr->count_shapes);
			tbr->ysh_min = (coord*)aligned_alloc(16, sizeof(coord) * tbr->count_shapes);
			tbr->ysh_max = (coord*)aligned_alloc(16, sizeof(coord) * tbr->count_shapes);
			indexer k = 0;
#ifdef OLD_LEAFS
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
#else
			// TO DO LEAFS
			indexer tj = tbr->leaf_number[0];
			tbr->offset[k] = 0;
			coord xsh_min, xsh_max, ysh_min, ysh_max;
			xsh_min = xsh_max = tbr->leaf_x[0];
			ysh_min = ysh_max = tbr->leaf_y[0];
			for (indexer j = 1; j < tbr->count_leafs; ++j) {
				if (tj != tbr->leaf_number[j]) {
					tbr->xsh_min[k] = xsh_min;
					tbr->xsh_max[k] = xsh_max;
					tbr->ysh_min[k] = ysh_min;
					tbr->ysh_max[k] = ysh_max;
					++k;
					tbr->offset[k] = j;
					xsh_min = xsh_max = tbr->leaf_x[j];
					ysh_min = ysh_max = tbr->leaf_y[j];
					tj = tbr->leaf_number[j];
					continue;
				}
				if (xsh_min > tbr->leaf_x[j])
					xsh_min = tbr->leaf_x[j];
				else if (xsh_max < tbr->leaf_x[j])
					xsh_max = tbr->leaf_x[j];
				if (ysh_min > tbr->leaf_y[j])
					ysh_min = tbr->leaf_y[j];
				else if (ysh_max < tbr->leaf_y[j])
					ysh_max = tbr->leaf_y[j];
			}
#endif // OLD_LEAFS
			// boundary of branch
			tbr->xsh_min[k] = xsh_min;
			tbr->xsh_max[k] = xsh_max;
			tbr->ysh_min[k] = ysh_min;
			tbr->ysh_max[k] = ysh_max;
			++k;
			//printf("SHAPES %u = %u (%u)\n", i, tbr->count_shapes, k);
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
		return m_nodes;
	}

	return NULL;
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

/// init empty root of tree
void init_root(const void* p) {
	m_nodes = (struct node*)aligned_alloc(16, sizeof(struct node));
	m_nodes->count_child_nodes = 1;
	alignas(16) struct branch *tbr = (struct branch*)aligned_alloc(16, sizeof(struct branch));
//#ifndef _WIN
	//*(m_nodes->child_node) = malloc(sizeof(struct branch*) * 1);
//#else
	m_nodes->child_node = (void**)aligned_alloc(16, sizeof(void*) * 1);
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
	m_nodes->center_child_node = (struct center_node_st*)aligned_alloc(16, sizeof(struct center_node_st) * 1);
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
#ifdef OLD_LEAFS
	tbr->leafs =  NULL;
#else
	// TO DO LEAFS
	tbr->leaf_x = NULL;
	tbr->leaf_y = NULL;
	tbr->leaf_number = NULL;
#endif // OLD_LEAFS
}

/// init empty root of tree v2
#ifdef OLD_LEAFS
void init_root2(struct node *nd, const void* p) {
#else
void init_root2(struct node *nd, coord x, coord y) {
#endif // OLD_LEAFS
	// nd = (struct node*)malloc(sizeof(struct node));
	nd->count_child_nodes = 1;
	alignas(16) struct branch *tbr = (struct branch*)aligned_alloc(16, sizeof(struct branch));
//#ifndef _WIN
	//*(nd->child_node) = malloc(sizeof(struct branch*) * 1);
//#else
	nd->child_node = (void**)aligned_alloc(16, sizeof(void*) * 1);
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
	nd->center_child_node = (struct center_node_st*)aligned_alloc(16, sizeof(struct center_node_st) * 1);
	nd->center_child_node[0].cy = 0.0;
	nd->center_child_node[0].cx = 0.0;
	nd->center_child_node[0].pos = NULL;

	// set first leaf
#ifdef OLD_LEAFS
	struct leaf *lf = (struct leaf*)p;

	nd->x1 = lf->x;
	nd->y1 = lf->y;
	nd->x2 = lf->x;
	nd->y2 = lf->y;
#else
	// TO DO LEAFS
	nd->x1 = x;
	nd->y1 = y;
	nd->x2 = x;
	nd->y2 = y;
#endif // OLD_LEAFS

	//struct branch *br = (struct branch*)m_nodes->child_node[0];
	tbr->count_leafs = 0;
	//br->count_merged_leafs = NULL;
	tbr->merge_next_leaf = NULL;
	tbr->curr_mem_pos = 0;
	tbr->alloc_mem_times = 0;
#ifdef OLD_LEAFS
	tbr->leafs = NULL;
#else
	// TO DO LEAFS
	tbr->leaf_x = NULL;
	tbr->leaf_y = NULL;
	tbr->leaf_number = NULL;
#endif // OLD_LEAFS
}

/// delete root
void del_root()
{
	//const size_t mem_size = 128;
	//unsigned mem_offset = 1;
	//struct node **stack_node = (struct node**)malloc(sizeof(struct node*) * mem_size * mem_offset);
	//unsigned stack_pos = 0;
	//indexer *stack_idx = (indexer*)malloc(sizeof(indexer) * mem_size * mem_offset);

	alignas(16) struct node *nd = m_nodes;

	if (nd->is_last_node) {
		for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
			alignas(16) struct branch *br = (struct branch*)(nd->child_node)[i];
			if (br->merge_next_leaf)
				_aligned_free(br->merge_next_leaf);
#ifdef OLD_LEADS
			if (br->leafs)
				_aligned_free(br->leafs);
#else
#endif // OLD_LEAFS
			if (br->length)
				_aligned_free(br->length);
		}
		_aligned_free((struct branch*)nd->child_node[0]);
		_aligned_free(nd->child_node);
		_aligned_free(nd);
		//free(m_nodes->child_node);
		//free(m_nodes);
	}
	else {
		//struct node* nd1 = NULL;
		indexer i = 0;
		alignas(16) struct branch *first_branch = NULL;
		//alignas(16) struct node* first_node = NULL;
		alignas(16) struct node* stack_first_node[64];
		//alignas(16) struct node* stack_child_node[64];
		for (unsigned i = 0; i < 64; ++i) {
			stack_first_node[i] = NULL;
			//stack_child_node[i] = NULL;
		}
		struct node *stack_node[64];
		int stack_pos = 0;
		indexer stack_idx[64];

		while (i < nd->count_child_nodes) {
			if (!stack_first_node[stack_pos] || nd < stack_first_node[stack_pos])
				stack_first_node[stack_pos] = nd;
			if (nd->is_last_node) {
				/*if (!first_node || nd < first_node)
					first_node = nd; */
				for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
					struct branch *br = (struct branch*)(nd->child_node[j]);
					if (!first_branch || br < first_branch)
						first_branch = br;
					_aligned_free(br->merge_next_leaf);
					//free(br->center);
#ifdef OLD_LEAFS
					_aligned_free(br->leafs);
#else
#endif // OLD_LEAFS
					_aligned_free(br->length);
					// boundary
					_aligned_free(br->xsh_max);
					_aligned_free(br->xsh_min);
					_aligned_free(br->ysh_max);
					_aligned_free(br->ysh_min);
					_aligned_free(br->offset);
					// leafs
					_aligned_free(br->leaf_x);
					_aligned_free(br->leaf_y);
					_aligned_free(br->leaf_number);
				}
				_aligned_free(nd->child_node);
				_aligned_free(nd->center_child_node);
				//free(nd);
				// return from stack
				while (stack_pos > 0) {
					stack_pos--;
					nd = stack_node[stack_pos];
					i = stack_idx[stack_pos] + 1;

					if (i < nd->count_child_nodes) {
						// insert to stack
						stack_idx[stack_pos] = i;
						stack_node[stack_pos] = nd;
						stack_pos++;
						nd = (struct node*)nd->child_node[i];
						i = 0;
						break;
					} else {
						//free(nd->child_node[0]);
						_aligned_free(nd->child_node);
						//free(nd->center_child_node);
						//free(nd);
					}
				}
			} else {
				// insert to stack
				stack_idx[stack_pos] = i;
				stack_node[stack_pos] = nd;
				stack_pos++;
				_aligned_free(nd->center_child_node);
				nd = (struct node*)nd->child_node[i];
			}
		}
//		return;
		_aligned_free(first_branch);
		//_aligned_free(first_node);
		for (unsigned k = 0; k < 64; ++k) {
			if (stack_first_node[k])
				_aligned_free(stack_first_node[k]);
			else
				break;
		}
	}
}

/// add leafs v2
bool add_leafs2(struct branch *br, struct leaf *lf, indexer count_of_leafs)
{
	// temp for count
	// unsigned tmp = 0; // br->count_leafs;
	// struct branch *br = m_nodes[i].child_node;
#ifdef OLD_LEAFS
	if (br->leafs == NULL) {
		br->leafs = (struct leaf*)aligned_alloc(16, sizeof(struct leaf) * m_max_added_leafs_in_iter);
		br->merge_next_leaf = (bool*)aligned_alloc(16, sizeof(bool) * m_max_added_leafs_in_iter);
		br->alloc_mem_times = 1;
	} else {
		if (br->count_leafs + count_of_leafs > br->alloc_mem_times * m_max_added_leafs_in_iter) {
			++(br->alloc_mem_times);
			br->leafs = (struct leaf*)_aligned_realloc(br->leafs, sizeof(struct leaf) * br->alloc_mem_times * m_max_added_leafs_in_iter, 16);
			br->merge_next_leaf = (bool*)_aligned_realloc(br->merge_next_leaf, sizeof(bool) * br->alloc_mem_times * m_max_added_leafs_in_iter, 16);
		}
	}
	//br->count_merged_leafs = (unsigned*)malloc(sizeof(unsigned) * 1);
	//br->merge_next_leaf = (bool*)malloc(sizeof(bool) * count_of_leafs);
	br->count_leafs += count_of_leafs;
	memcpy(br->leafs + br->curr_mem_pos, lf, sizeof(struct leaf) * count_of_leafs);
#else
	// TO DO LEAFS
	if (br->leaf_x == NULL) {
		br->leaf_x = (coord*)aligned_alloc(16, sizeof(coord) * m_max_added_leafs_in_iter);
		br->leaf_y = (coord*)aligned_alloc(16, sizeof(coord) * m_max_added_leafs_in_iter);
		br->leaf_number = (indexer*)aligned_alloc(16, sizeof(indexer) * m_max_added_leafs_in_iter);
		br->merge_next_leaf = (bool*)aligned_alloc(16, sizeof(bool) * m_max_added_leafs_in_iter);
		br->alloc_mem_times = 1;
	}
	else {
		if (br->count_leafs + count_of_leafs > br->alloc_mem_times * m_max_added_leafs_in_iter) {
			++(br->alloc_mem_times);
			br->leaf_x = (coord*)_aligned_realloc(br->leaf_x, sizeof(coord) * br->alloc_mem_times * m_max_added_leafs_in_iter, 16);
			br->leaf_y = (coord*)_aligned_realloc(br->leaf_y, sizeof(coord) * br->alloc_mem_times * m_max_added_leafs_in_iter, 16);
			br->leaf_number = (indexer*)_aligned_realloc(br->leaf_number, sizeof(indexer) * br->alloc_mem_times * m_max_added_leafs_in_iter, 16);
			br->merge_next_leaf = (bool*)_aligned_realloc(br->merge_next_leaf, sizeof(bool) * br->alloc_mem_times * m_max_added_leafs_in_iter, 16);
		}
	}

	br->count_leafs += count_of_leafs;
	// memcpy(br->leafs + br->curr_mem_pos, lf, sizeof(struct leaf) * count_of_leafs);
	for (indexer i = 0; i < count_of_leafs; ++i) {
		br->leaf_x[br->curr_mem_pos + i] = lf[i].x;
		br->leaf_y[br->curr_mem_pos + i] = lf[i].y;
		br->leaf_number[br->curr_mem_pos + i] = lf[i].number;
	}
#endif // OLD_LEAFS
	unsigned i = br->curr_mem_pos;
	for (; i < br->count_leafs - 1; ++i) {
		(br->merge_next_leaf)[i] = true; // count_of_leafs
	}
	(br->merge_next_leaf)[i] = false;
	br->curr_mem_pos = br->count_leafs;
	return true;
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
	alignas(16) struct branch *br = (struct branch*)(nd->child_node[0]);
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
		alignas(16) struct boudary *positions1 = (struct boudary*)aligned_alloc(16, sizeof(struct boudary) * size_separate);
		alignas(16) struct boudary *positions2 = (struct boudary*)aligned_alloc(16, sizeof(struct boudary) * size_separate);
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
		alignas(16) struct branch *br1 = (struct branch*)aligned_alloc(16, sizeof(struct branch) * size_separate);
		indexer total_count_shapes = 0;
		for (unsigned i = 0; i < size_separate; ++i) {
			pos_start = positions1[i].idx1;
			//center_start = &(center_start[pos_start]);
			cshapes = positions1[i].idx2 - positions1[i].idx1 + 1;
			tcleafs = positions1[i].count_leafs;

#ifdef OLD_LEAFS
			br1[i].leafs = (struct leaf*)aligned_alloc(16, sizeof(struct leaf) * tcleafs);
#else
			// TO DO LEAFS
			br1[i].leaf_x = (coord*)aligned_alloc(16, sizeof(coord) * tcleafs);
			br1[i].leaf_y = (coord*)aligned_alloc(16, sizeof(coord) * tcleafs);
			br1[i].leaf_number = (indexer*)aligned_alloc(16, sizeof(indexer) * tcleafs);
#endif // OLD_LEAFS
			br1[i].merge_next_leaf = (bool*)aligned_alloc(16, sizeof(bool) * tcleafs);

			br1[i].alloc_mem_times = (unsigned)-1;
			br1[i].curr_mem_pos = (unsigned)-1;
			br1[i].count_shapes = cshapes;
			br1[i].count_leafs = tcleafs;
			indexer count_of_leafs = 0;
			for (indexer j = 0; j < cshapes; ++j) {
#ifdef OLD_LEAFS
				memcpy(br1[i].leafs + count_of_leafs, br->center[j + total_count_shapes].pos_leaf, sizeof(struct leaf) * br->center[j + total_count_shapes].count_leafs);
#else
				// TO DO LEAFS
				memcpy(br1[i].leaf_x + count_of_leafs, (void*)&(br->leaf_x[br->center[j + total_count_shapes].pos_leaf]), sizeof(coord) * br->center[j + total_count_shapes].count_leafs);
				memcpy(br1[i].leaf_y + count_of_leafs, (void*)&(br->leaf_y[br->center[j + total_count_shapes].pos_leaf]), sizeof(coord) * br->center[j + total_count_shapes].count_leafs);
				memcpy(br1[i].leaf_number + count_of_leafs, (void*)&(br->leaf_number[br->center[j + total_count_shapes].pos_leaf]), sizeof(indexer) * br->center[j + total_count_shapes].count_leafs);
#endif // OLD_LEAFS
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
#ifdef OLD_LEAFS
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
#else
			// TO DO LEAFS
			coord xmi = br1[i].leaf_x[0];
			coord xma = xmi;
			coord ymi = br1[i].leaf_y[0];
			coord yma = ymi;
			for (indexer j = 1; j < br1[i].count_leafs; ++j) {
				if (br1[i].leaf_x[j] < xmi)
					xmi = br1[i].leaf_x[j];
				if (br1[i].leaf_x[j] > xma)
					xma = br1[i].leaf_x[j];
				if (br1[i].leaf_y[j] < ymi)
					ymi = br1[i].leaf_y[j];
				if (br1[i].leaf_y[j] > yma)
					yma = br1[i].leaf_y[j];
			}
#endif // OLD_LEAFS
			br1[i].x_min = xmi;
			br1[i].x_max = xma;
			br1[i].y_min = ymi;
			br1[i].y_max = yma;
		}

		// delete old branch and exchange on curernt branches
		_aligned_free(br->center);
#ifdef OLD_LEAFS
		_aligned_free(br->leafs);
#else
		// TO DO LAEFS
		_aligned_free(br->leaf_x);
		_aligned_free(br->leaf_y);
		_aligned_free(br->leaf_number);
#endif // OLD_LEAFS
		_aligned_free(br->merge_next_leaf);
		_aligned_free(br);
		_aligned_free(nd->child_node);
		nd->child_node = (void**)aligned_alloc(16, sizeof(void*) * size_separate);
		_aligned_free(nd->center_child_node); // ???
		nd->center_child_node = (struct center_node_st*)aligned_alloc(16, sizeof(struct center_node_st) * size_separate);
		for (indexer i = 0; i < size_separate; ++i) {
			nd->child_node[i] = &(br1[i]);
			// centers
			nd->center_child_node[i].cx = (br1[i].x_max - br1[i].x_min) / (coord)2.0 + br1[i].x_min;
			nd->center_child_node[i].cy = (br1[i].y_max - br1[i].y_min) / (coord)2.0 + br1[i].y_min;
			nd->center_child_node[i].pos = &(br1[i]);
		}
		nd->count_child_nodes = (unsigned)size_separate;

		_aligned_free(positions1);
		_aligned_free(positions2);

	//	return nd;
	}

	return nd;
}

/// create threads for first node
bool create_first_thread(struct node* nd, struct leaf* leafs, unsigned *offsets_leafs, unsigned count_cpus, unsigned count_leafs)
{
	// http://learnc.info/c/pthreads_create_and_join.html
	// https://msdn.microsoft.com/en-us/library/kdzttdcb.aspx
	// pthread_create for linux, _beginthread for windows

#ifndef _WIN
	// handler thread
	alignas(16) pthread_t *ptr1 = (pthread_t*)aligned_alloc(16, sizeof(pthread_t) * count_cpus);
	//unsigned tun[32];
	alignas(16) struct first_thr_st *ftst = (struct first_thr_st*)aligned_alloc(16, sizeof(struct first_thr_st) * count_cpus);
	alignas(16) unsigned *idx = (unsigned*)aligned_alloc(16, sizeof(unsigned) * (count_cpus + 1));

	// prepare for separate
	unsigned offset = (unsigned)ceil((double)count_leafs / count_cpus);
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

	_aligned_free(idx);
	_aligned_free(ftst);
	_aligned_free(ptr1);

#else
	// handler thread
	alignas(16) uintptr_t *ptr1 = (uintptr_t*)aligned_alloc(16, sizeof(uintptr_t) * count_cpus);
	//unsigned tun[32];
	alignas(16) struct first_thr_st *ftst = (struct first_thr_st*)aligned_alloc(16, sizeof(struct first_thr_st) * count_cpus);
	alignas(16) unsigned *idx = (unsigned*)aligned_alloc(16, sizeof(unsigned) * (count_cpus + 1));

	// prepare for separate
	unsigned offset = (unsigned)ceil((double)count_leafs / count_cpus);
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

	_aligned_free(idx);
	_aligned_free(ftst);
	_aligned_free(ptr1);
#endif // _WIN
	return false;
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
#ifdef OLD_LEAFS
	init_root2(ftst->node_, &(ftst->leafs_)[ftst->start_pos_leafs]);
#else
	init_root2(ftst->node_, (ftst->leafs_)[ftst->start_pos_leafs].x, (ftst->leafs_)[ftst->start_pos_leafs].y);
#endif // OLD_LEAFS

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

	br->center = (struct center_st2*)aligned_alloc(16, sizeof(struct center_st2) * count_shapes);
	// hope that this is the fitst time and alloc memory
	//////////////////br->count_merged_pos_leafs = (unsigned*)malloc(sizeof(unsigned) * count_shapes);
	// br->center.pos_leaf = (struct leaf**)malloc(sizeof(struct leaf*) * count_shapes);
	//br->center->cx = (coord*)malloc(sizeof(coord) * count_shapes);
	//br->center->cy = (coord*)malloc(sizeof(coord) * count_shapes);

	unsigned j = 0; // count leafs in union
	unsigned k = 0; // index of shape
#ifdef OLD_LEAFS
	br->center[k].pos_leaf = &(br->leafs[0]);
#else
	// TO DO LEAFS
	br->center[k].pos_leaf = 0;
#endif // OLD_LEAFS
	//j++;
	for (unsigned i = 0; i < br->count_leafs; ++i) {
		if (!br->merge_next_leaf[i]) {
			// calculate center
			tx = 0.0;
			ty = 0.0;
			for (unsigned t1 = i - j; t1 <= i; ++t1) {
#ifdef OLD_LEAFS
				tx += br->leafs[t1].x;
				ty += br->leafs[t1].y;
#else
				// TO DO LEAFS
				tx += br->leaf_x[t1];
				ty += br->leaf_y[t1];
#endif // OLD_LEAFS
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
#ifdef OLD_LEAFS
			br->center[k].pos_leaf = &(br->leafs[i + 1]);
#else
			// TO DO LEAFS
			br->center[k].pos_leaf = i + 1;
#endif // OLD_LEAFS
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
		return nd;
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
	alignas(16) struct boudary *positions1 = (struct boudary*)aligned_alloc(16, sizeof(struct boudary) * size_separate);
	alignas(16) struct boudary *positions2 = (struct boudary*)aligned_alloc(16, sizeof(struct boudary) * size_separate);
	//unsigned tcbranches = cbranches;
	unsigned brinscope = cbranches;

	positions2[0].idx1 = pos_start;
	positions2[0].idx2 = cbranches - 1;
	//indexer j;
	indexer k;
	//indexer ccount;

	bool pos_pos = false; // false - used positions1, true = used positions2

	for (unsigned i = 2; i <= size_separate; i = i << 2) {
		//tcbranches = cbranches / i; // ???
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
		//tcbranches = cbranches / (i * 2);  // ???
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
	alignas(16) struct branch* tbr1 = (struct branch*)aligned_alloc(16, sizeof(struct branch) * nd->count_child_nodes);
	for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
		memcpy(&(tbr1[i]), nd->center_child_node[i].pos, sizeof(struct branch));
	}

	// prepare new nodes
	alignas(16) struct node *nd1 = (struct node*)aligned_alloc(16, sizeof(struct node) * size_separate);
	//__declspec(align(16)) void** pointer_void = (void**)aligned_alloc(16, sizeof(void*) * size_separate);
	//indexer total_count_br = 0;
	for (unsigned i = 0; i < size_separate; ++i) {
#ifdef OLD_LEAFS
		init_root2(&(nd1[i]), &(((struct branch*)(nd->child_node[0]))->leafs[0]));
#else
		// TO DO LEAFS
		//init_root2(&(nd1[i]), ((struct branch*)(nd->child_node[0]))->leaf_x[0], ((struct branch*)(nd->child_node[0]))->leaf_y[0]);
		init_root2(&(nd1[i]), tbr1[positions1[i].idx1].x_min, tbr1[positions1[i].idx1].y_min);
#endif // OLD_LEAFS
		// assign branches to node
		_aligned_free((struct branch*)nd1[i].child_node[0]);
		_aligned_free(nd1[i].child_node);
		nd1[i].child_node = (void**)aligned_alloc(16, sizeof(void*) * (positions1[i].idx2 - positions1[i].idx1 + 1));
		//nd1[i].child_node = &(pointer_void[positions1[i].idx1]);
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

	// free old branch, i.e. exchanged to tbr1
	_aligned_free(nd->child_node[0]);
	// free branches and nodes
	_aligned_free(nd->child_node);
	nd->child_node = (void**)aligned_alloc(16, sizeof(void*) * size_separate);
	// free centers
	_aligned_free(nd->center_child_node);
	nd->center_child_node = (struct center_node_st*)aligned_alloc(16, sizeof(struct center_node_st) * size_separate);
	for (indexer i = 0; i < size_separate; ++i) {
		nd->child_node[i] = &(nd1[i]);
		// centers
		nd->center_child_node[i].cx = (nd1[i].x2 - nd1[i].x1) / (coord)2.0 + nd1[i].x1;
		nd->center_child_node[i].cy = (nd1[i].y2 - nd1[i].y1) / (coord)2.0 + nd1[i].y1;
		nd->center_child_node[i].pos = &(nd1[i]);
	}
	nd->count_child_nodes = (unsigned)size_separate;
	nd->is_last_node = false;

	_aligned_free(positions1);
	_aligned_free(positions2);

	return nd;
////////////////////////////////////////////////////////////////

	return NULL;
}

/// separate nodes
struct node* separate_nodes(struct node* nd)
{
	if (nd->is_last_node) {
		return nd;
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
		alignas(16) struct boudary *positions1 = (struct boudary*)aligned_alloc(16, sizeof(struct boudary) * size_separate);
		alignas(16) struct boudary *positions2 = (struct boudary*)aligned_alloc(16, sizeof(struct boudary) * size_separate);
		//unsigned tcnodes = cnodes;
		unsigned ndinscope = cnodes;

		positions2[0].idx1 = pos_start;
		positions2[0].idx2 = cnodes - 1;
		//indexer j;
		indexer k;
		//indexer ccount;

		bool pos_pos = false; // false - used positions1, true = used positions2

		for (unsigned i = 2; i <= size_separate; i = i << 2) {
			//tcnodes = cnodes / i; // ???
			pos_pos = false;
			for (k = 0; k < i >> 1; ++k) {
				pos_start = positions2[k].idx1;
				center_start = &(nd->center_child_node[pos_start]);
				ndinscope = positions2[k].idx2 - positions2[k].idx1 + 1;
				//ccount = 0;

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
			//tcnodes = cnodes / (i * 2);  // ???
			for (k = 0; k < i; ++k) {
				pos_start = positions1[k].idx1;
				center_start = &(nd->center_child_node[pos_start]);
				ndinscope = positions1[k].idx2 - positions1[k].idx1 + 1;
				//ccount = 0;
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
		alignas(16) struct node* tnd1 = (struct node*)aligned_alloc(16, sizeof(struct node) * nd->count_child_nodes);
		for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
			memcpy(&(tnd1[i]), nd->center_child_node[i].pos, sizeof(struct node));
		}

		// prepare new nodes
		alignas(16) struct node *nd1 = (struct node*)aligned_alloc(16, sizeof(struct node) * size_separate);
		//indexer total_count_br = 0;
		for (unsigned i = 0; i < size_separate; ++i) {
#ifdef OLD_LEAFS
			struct leaf lf;
			lf.x = 0.0;
			lf.y = 0.0;
			lf.number = 0;
			init_root2(&(nd1[i]), &lf);
#else
			// TO DO LEAFS
			init_root2(&(nd1[i]), 0.0, 0.0);
#endif
			// assign branches to node
			_aligned_free(/*(struct branch*)*/nd1[i].child_node[0]);
			_aligned_free(nd1[i].child_node);
			nd1[i].child_node = (void**)aligned_alloc(16, sizeof(void*) * (positions1[i].idx2 - positions1[i].idx1 + 1));
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
		_aligned_free(nd->child_node[0]);
		_aligned_free(nd->child_node);
		nd->child_node = (void**)aligned_alloc(16, sizeof(void*) * size_separate);
		// free centers
		_aligned_free(nd->center_child_node);
		nd->center_child_node = (struct center_node_st*)aligned_alloc(16, sizeof(struct center_node_st) * size_separate);
		for (indexer i = 0; i < size_separate; ++i) {
			nd->child_node[i] = &(nd1[i]);
			// centers
			nd->center_child_node[i].cx = (nd1[i].x2 - nd1[i].x1) / (coord)2.0 + nd1[i].x1;
			nd->center_child_node[i].cy = (nd1[i].y2 - nd1[i].y1) / (coord)2.0 + nd1[i].y1;
			nd->center_child_node[i].pos = &(nd1[i]);
		}
		nd->count_child_nodes = (unsigned)size_separate;
		nd->is_last_node = false;

		_aligned_free(positions1);
		_aligned_free(positions2);
	}

	// update boundary of top node
	coord x_max, y_max, x_min, y_min, tx_max, ty_max, tx_min, ty_min;
	if (nd->is_last_node) {
		// branch
		x_min = ((struct branch*)nd->child_node[0])->x_min;
		x_max = ((struct branch*)nd->child_node[0])->x_max;
		y_min = ((struct branch*)nd->child_node[0])->y_min;
		y_max = ((struct branch*)nd->child_node[0])->y_max;
	}
	else {
		// node
		x_min = ((struct node*)nd->child_node[0])->x1;
		x_max = ((struct node*)nd->child_node[0])->x2;
		y_min = ((struct node*)nd->child_node[0])->y1;
		y_max = ((struct node*)nd->child_node[0])->y2;
	}
	for (indexer i = 1; i < nd->count_child_nodes; ++i) {
		if (nd->is_last_node) {
			// branch
			tx_min = ((struct branch*)nd->child_node[i])->x_min;
			tx_max = ((struct branch*)nd->child_node[i])->x_max;
			ty_min = ((struct branch*)nd->child_node[i])->y_min;
			ty_max = ((struct branch*)nd->child_node[i])->y_max;
		}
		else {
			// node
			tx_min = ((struct node*)nd->child_node[i])->x1;
			tx_max = ((struct node*)nd->child_node[i])->x2;
			ty_min = ((struct node*)nd->child_node[i])->y1;
			ty_max = ((struct node*)nd->child_node[i])->y2;
		}
		if (tx_min < x_min)
			x_min = tx_min;
		if (tx_max > x_max)
			x_max = tx_max;
		if (ty_min < y_min)
			y_min = ty_min;
		if (ty_max > y_max)
			y_max = ty_max;
	}
	nd->x1 = x_min;
	nd->x2 = x_max;
	nd->y1 = y_min;
	nd->y2 = y_max;

	return nd;
	////////////////////////////////////////////////////////////////

	return NULL;
}
