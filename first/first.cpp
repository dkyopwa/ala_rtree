#include "pch.h"
#include <stdio.h>
#include <malloc.h>
#ifndef _WIN
	#include <sys/time.h>
#else
	#include "time.h"
	#include <conio.h>
#endif
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "sort.h"
#include "log.h"


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

#define indexer	unsigned

const double PI = 3.1415926535897932384626433832795;
const unsigned MAX_ITEMS_IN_NODE = 29950;
const unsigned MAX_NODES = 10;

struct leaf {
	coord x;
	coord y;
	indexer number;	// number of item in the global scope
#ifndef _WIN
} __attribute__ ((aligned (16)));
#else
};
#endif

struct branch {
	unsigned count_leafs;
	struct leaf* leafs;
	//unsigned *count_merged_leafs; // 0 - not leaf, 1 and more - count of point in the branch
	bool *merge_next_leaf;
#ifndef _WIN
} __attribute__ ((aligned (16)));
#else
};
#endif

struct node {
	// boundary
	coord x1;
	coord y1;
	coord x2;
	coord y2;
	void** child_node; // may be node or branch
	bool is_last_node;
	unsigned count_child_nodes;
#ifndef _WIN
} __attribute__ ((aligned (16)));
#else
};
#endif

// globale variables
struct node* m_nodes = NULL;

/// generate data for test
struct leaf* generate(unsigned *count, unsigned **offsets_leafs);
void lprintf(const char *text);
void init_root(const void* p);
void del_root();
void add_leafs(struct branch *br, struct leaf *lf, unsigned count_of_leafs);
bool add(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs);
bool add_leaf_in_boundary(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs);
bool add_leaf_out_boundary(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs);
struct node* check_separate(struct node* nd);
struct node* separate(struct node* nd);
struct node* separate_leafs(struct node* nd, unsigned count_senters_items, unsigned *pos_idx); //_x, unsigned *pos_idx_y, bool separate_by_x);

/// main function
int main()
{
	lprintf("Hello");

	unsigned count_of_leafs = 0; // old: must to devide 3 without
	unsigned *offsets_leafs = NULL;
	struct leaf* lll = generate(&count_of_leafs, &offsets_leafs);

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

	char tch[1024] = {0};
///////////////	sprintf(tch, "count = %d, x_min = %f, x_max = %f, y_min = %f, y_max = %f", m_nodes->count_child_nodes, m_nodes->x1, m_nodes->x2, m_nodes->y1, m_nodes->y2);
	lprintf(tch);

	if (lll)
		free(lll);
	//if (m_nodes)
	//	free(m_nodes);
	del_root();

	lprintf("Done");
#ifdef _WIN
	_getch();
#endif
	return 0;
}



/// generate test data
struct leaf* generate(unsigned *count, unsigned **offsets_leafs)
{
	FILE *f1;
	errno_t t = fopen_s(&f1, "c:/projects/tmp/1/1.bin", "rb");

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

	struct leaf* res = NULL;
	res = (struct leaf*)malloc(sizeof(struct leaf) * count_leafs);
	if (!res)
		return NULL;

	// struct data *dts = (struct data*)malloc(sizeof(struct data) * count1);
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
		//dts[i].pts = (struct point*)malloc(sizeof(struct point) * cnt);
		//st = fread(dts[i].pts, sizeof(double), cnt * 2, f1);
	}

	fclose(f1);

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
#ifndef _WIN
	//*(m_nodes->child_node) = malloc(sizeof(struct branch*) * 1);
#else
	m_nodes->child_node = (void**)malloc(sizeof(void*) * 1);
#endif
	(m_nodes->child_node)[0] = tbr;
	/*m_nodes->is_node = (bool*)malloc(sizeof(bool) * 1);
	m_nodes->is_node[0] = false;
	*/
	m_nodes->is_last_node = true;

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
	tbr->leafs =  NULL;
}

void del_root()
{
	struct branch *br = (struct branch*)m_nodes->child_node;
	if (br->merge_next_leaf)
		free(br->merge_next_leaf);
	if (br->leafs)
		free(br->leafs);
	//free(m_nodes->is_node);
	free(m_nodes->child_node);
	free(m_nodes);
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

/// add item to tree
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
		add_leafs(&((struct branch *)(nd->child_node))[idx_node], l, count_of_leafs);
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

bool add_leaf_out_boundary(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs)
{
	// TO DO +++
	// not in the boundary
	unsigned idx_node = 0;
	bool flag = false;
	// is leaf
	struct leaf *l = (struct leaf*)p;

	// need to expand boundary
	coord min_square = 999999.0;
	coord min_fake_x1 = 999999.0, min_fake_y1 = 999999.0, min_fake_x2 = 999999.0, min_fake_y2 = 999999.0;
	// select node to expand
	//unsigned idx = 0;
	coord fake_x1, fake_y1, fake_x2, fake_y2;
	for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
		struct node* tnd = &((struct node*)(nd->child_node))[j];
		for (unsigned k = 0; k < count_of_leafs; ++k) {
			if (l[k].x < tnd->x1 && l[k].x < tnd->x2) {
				// expand left boundary
				fake_x1 = l[k].x;
				if (!k)
					fake_x2 = tnd->x2;
			} else if (l[k].x > tnd->x1 && l[k].x > tnd->x2) {
				// expand right boundary
				if (!k)
					fake_x1 = tnd->x1;
				fake_x2 = l[k].x;
			} else if (!k) {
				fake_x1 = tnd->x1;
				fake_x2 = tnd->x2;
			}
			if (l[k].y < tnd->y1 && l[k].y < tnd->y2) {
				// expand left boundary
				fake_y1 = l[k].y;
				if (!k)
					fake_y2 = tnd->y2;
			} else if (l[k].y > tnd->y1 && l[k].y > tnd->y2) {
				// expand right boundary
				if (!k)
					fake_y1 = tnd->y1;
				fake_y2 = l[k].y;
			} else if (!k) {
				fake_y1 = tnd->y1;
				fake_y2 = tnd->y2;
			}
		}
		coord square = fabs(fake_x1 - fake_x2) * fabs(fake_y1 - fake_y2);
		char tch[1024] = {0};
///////////		sprintf(tch, "square = %f, x1 = %f, x2 = %f, y1 = %f, y2 = %f, count_of_leafs = %d", square, fake_x1, fake_x2, fake_y1, fake_y2, count_of_leafs);
///////////		lprintf(tch);
		if (square < min_square) {
			idx_node = j;
			min_square = square;
			min_fake_x1 = fake_x1;
			min_fake_y1 = fake_y1;
			min_fake_x2 = fake_x2;
			min_fake_y2 = fake_y2;
		}
	}

	// recusive to childen nodes
	if (!nd->is_last_node) {
		//for  (unsigned i = 0; i < nd->count_nodes; ++i) {
			if (add_leaf_out_boundary(&((struct node*)(nd->child_node))[idx_node], is_leaf, p, count_of_leafs)) {
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
		m_nodes[idx_node].x1 = min_fake_x1;
		m_nodes[idx_node].y1 = min_fake_y1;
		m_nodes[idx_node].x2 = min_fake_x2;
		m_nodes[idx_node].y2 = min_fake_y2;
		// add leaf
		if (nd->is_last_node)
			add_leafs(&((struct branch *)(nd->child_node))[idx_node], l, count_of_leafs);
		return true;
	}
	return false;
}

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
			q_sort_mas(xc, count_centers, pos_idx_x);
			/* f1 = fopen("2.txt", "w");
			if (f1) {
				for (unsigned i = 0; i < count_centers; ++i)
					fprintf(f1, "%d: xc = %f, yc = %f\n", pos_idx_x[i], xc[i], yc[i]);
				fclose(f1);
			} */
			q_sort_mas(yc, count_centers, pos_idx_y);
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
			min_max_sort(&area_x, &area_y);
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
				min_max_sort(&x1, &y1);
				min_max_sort(&x2, &y2);
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
			struct node *tnd = (struct node*)nd->child_node[i];
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

struct node* separate(struct node* nd)
{

	return NULL;
}