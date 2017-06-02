/*
	Project present like RTree
	Created by Vladimir Nedved 2017
	Apache License 2.0
*/
#pragma once

#ifndef FIRST_H_HEADERS
#define FIRST_H_HEADERS

typedef float	coord;
typedef unsigned	indexer;

typedef bool (*ret_callback)(indexer, void*);

const double PI = 3.1415926535897932384626433832795;
const unsigned MAX_ITEMS_IN_NODE = 1000;
const unsigned MAX_NODES = 4;
const unsigned MAX_ADDED_LEAFS_IN_ITER = 1000000; // for malloc and realloc memory

struct leaf {
	coord x;
	coord y;
	indexer number;	// number of item in the global scope
#ifndef _WIN
} __attribute__((aligned(16)));
#else
};
#endif

struct center_st {
	struct leaf** pos_leaf;
	// unsigned *count_merged_pos_leafs;
	coord *cx;
	coord *cy;
	unsigned count_shapes;
#ifndef _WIN
} __attribute__((aligned(16)));
#else
};
#endif

struct center_st2 {
	struct leaf* pos_leaf;
	//unsigned *count_merged_pos_leafs;
	coord cx;
	coord cy;
	indexer count_leafs;
	//unsigned count_shapes;
	// bool moved;
#ifndef _WIN
} __attribute__((aligned(16)));
#else
};
#endif

struct center_node_st {
	coord cx;
	coord cy;
	void *pos;
#ifndef _WIN
} __attribute__((aligned(16)));
#else
};
#endif

struct branch {
	unsigned count_leafs;
	unsigned curr_mem_pos;
	unsigned alloc_mem_times;
	struct leaf* leafs;
	//unsigned *count_merged_leafs; // 0 - not leaf, 1 and more - count of point in the branch
	// for find center
	// struct leaf** pos_leaf;
	//////////// unsigned *count_merged_pos_leafs;
	unsigned count_shapes;
	/* coord *cx;
	coord *cy;
	*/
	struct center_st2 *center;
	// end for find center
	bool *merge_next_leaf;
	// boundary, calculate after separate
	coord x_min;
	coord y_min;
	coord x_max;
	coord y_max;
	//struct center_node_st branch_center;
#ifndef _WIN
} __attribute__((aligned(16)));
#else
};
#endif

struct node {
	// boundary
	coord x1;
	coord y1;
	coord x2;
	coord y2;
	void** child_node; // may be node (is_last_node = false) or branch (is_last_node = true)
	unsigned count_child_nodes;
	bool is_last_node;

	struct center_node_st *center_child_node;
#ifndef _WIN
} __attribute__((aligned(16)));
#else
};
#endif

struct first_thr_st {
	struct node* node_;
	struct leaf* leafs_;
	unsigned *offsets_leafs_;
	unsigned count_leaf;
	unsigned start_pos_leafs;
#ifndef _WIN
} __attribute__((aligned(16)));
#else
};
#endif

struct boudary {
	indexer idx1;
	indexer idx2;
	indexer count_leafs;
#ifndef _WIN
} __attribute__((aligned(16)));
#else
};
#endif

inline float very_fast_sqrt(float f)
{
/*	__asm movss   xmm0, [f];
	__asm rsqrtss xmm0, xmm0;
	__asm rcpss   xmm0, xmm0;
	__asm movss[f], xmm0;
	*/
	return f;
}

/// generate data for test
struct leaf* generate(unsigned *count, unsigned **offsets_leafs, unsigned *count_shapes);
/// print in stdout text with time with milliseconds
void lprintf(const char *text);
/// initialization root
void init_root(const void* p);
/// initialization root v2
void init_root2(struct node *nd, const void* p);
/// deleting root
void del_root();
void add_leafs(struct branch *br, struct leaf *lf, unsigned count_of_leafs);
bool add(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs);
bool add_leaf_in_boundary(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs);
bool add_leaf_out_boundary(struct node* nd, bool is_leaf, const void* p, unsigned count_of_leafs);
struct node* check_separate(struct node* nd);
struct node* separate(struct node* nd);
struct node* separate_branches(struct node* nd);
struct node* separate_nodes(struct node* nd);
//struct node* separate2(struct node* nd);
struct node* separate_leafs(struct node* nd, unsigned count_senters_items, unsigned *pos_idx); //_x, unsigned *pos_idx_y, bool separate_by_x);
bool create_first_thread(struct node* nd, struct leaf* leafs, unsigned *offsets_leafs, unsigned count_cpus, unsigned count_leafs);
void first_thread(void *params);
#ifndef _WIN
void* first_thread_v2(void *params);
#else
void first_thread_v2(void *params);
#endif
void print_file_svg(struct node *nds, unsigned count_nodes, unsigned count_shapes, const char *file_name);
bool find_centers(const struct node *nd, const unsigned count_shapes);
bool find_branch_centers(struct node *nd);
int x_cmp(const void* a, const void* b);
int y_cmp(const void* a, const void* b);
indexer search_rect(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, ret_callback callback, void *arg);
indexer* search_rects(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max);
indexer search_circle(struct node *nd, coord radius, ret_callback callback, void *arg);
indexer* search_circles(struct node *nd, coord radius);
indexer search_point(struct node *nd, coord x, coord y, coord radius);
void try_find(struct leaf* lll, unsigned count_of_leafs);


#endif //FIRST_H_HEADERS
