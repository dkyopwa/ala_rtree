#include <stdio.h>
//#include <malloc.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include "unimem.h"
#include "first.h"

#ifndef _WIN
	#define __int64 long long
#endif //_WIN

//#define MINIMAL_DEBUG2

struct point {
	coord x;
	coord y;
	point(coord xx, coord yy) : x(xx), y(yy) {};
	point() {};
};

struct calc_data {
	struct branch *br;
	indexer idx;
	union {
		coord radius;
		coord dist;
	};
	struct point center;
#ifdef CALC_POINT
	// for search point
	indexer curr_idx;
#endif // CALC POINT

	calc_data(): center(0, 0) {};
};

/*
	calculate distance from point to line
*/
coord distance2(struct point *p, struct point *line_p0, struct point *line_p1)
{
	coord vx, vy, wx, wy, c1, c2, b, pbx, pby;

	vx = line_p1->x - line_p0->x;
	vy = line_p1->y - line_p0->y;
	wx = p->x - line_p0->x;
	wy = p->y - line_p0->y;
	c1 = vx * wx + vy * wy;

	if (c1 <= 0) {
		//coord t1 = p->x - line_p0->x;
		//coord t2 = p->y - line_p0->y;
		coord t1 = wx;
		coord t2 = wy;
		return sqrt(t1 * t1 + t2 * t2);
	}

	c2 = vx * vx + vy * vy;
	if (c2 <= c1) {
		//return sqrt(pow(fabs(p->x - line_p1->x), 2) + pow(fabs(p->y - line_p1->y), 2));
		coord t1 = p->x - line_p1->x;
		coord t2 = p->y - line_p1->y;
		return sqrt(t1 * t1 + t2 * t2);
	}

	b = c1 / c2;
	pbx = line_p0->x + b * vx;
	pby = line_p0->y + b * vy;

	//return sqrt(pow(fabs(p->x - pbx), 2) + pow(fabs(p->y - pby), 2));
	coord t1 = p->x - pbx;
	coord t2 = p->y - pby;
	return sqrt(t1 * t1 + t2 * t2);
}

/*
	adding number of items from branch
*/
void add_branch(/*in*/struct branch* br, /*in/out*/size_t* mem_size, /*in/out*/size_t* count_mem, /*in/out*/indexer* idx, /*in/out*/indexer** idxs)
{
	// check enought memory
	indexer tmp = *idx + br->count_shapes;
	if (tmp > *mem_size * *count_mem) {
		__int64 tmp1 = *mem_size * *count_mem - tmp;
		if (tmp1 < 0) {
			*count_mem += (indexer)ceil(tmp1 * (-1.0) / *mem_size);
			*idxs = (indexer*)_aligned_realloc(*idxs, sizeof(indexer) * *mem_size * *count_mem, 16);
		}
	}

	for (indexer i1 = 0; i1 < br->count_shapes; ++i1) {
		// check memory and store index
		/*if (*idx > *mem_size * *count_mem) {
			(*count_mem)++;
			*idxs = (indexer*)_aligned_realloc(*idxs, sizeof(indexer) * *mem_size * *count_mem, 16);
		} */
		// store index of current sergment from shape
		(*idxs)[*idx] = br->leaf_number[br->offset[i1]];
		(*idx)++;

		indexer to_end;
		if (i1 == br->count_shapes - 1)
			to_end = br->count_leafs;
		else
			to_end = br->offset[i1 + 1];
	}
}

/*
	adding numbers of items from node
*/
void add_nodes(/*in*/struct node* nd, /*in/out*/size_t* mem_size, /*in/out*/size_t* count_mem, /*in/out*/indexer* idx, /*in/out*/indexer** idxs)
{
	indexer i = 0;
	struct node *stack_node[64];
	int stack_pos = 0;
	indexer stack_idx[64];

	while (i < nd->count_child_nodes) {
		if (nd->is_last_node) {
			for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
				struct branch *br = (struct branch*)(nd->child_node)[j];
				add_branch(br, mem_size, count_mem, idx, idxs);
			}
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
			}
		}
		else {
			stack_idx[stack_pos] = i;
			stack_node[stack_pos] = nd;
			stack_pos++;
			i = 0;
			nd = (struct node*)nd->child_node[i];
		}
	}
}

bool check_intersection(coord p1x, coord p1y, coord p2x, coord p2y, coord p3x, coord p3y, coord p4x, coord p4y)
{
	coord t0x = p2x - p1x;
	coord t0y = p2y - p1y;
	coord t1 = (p3x - p1x) * t0x + (p3y - p1y) * t0y;
	coord t2 = (p4x - p1x) * t0x + (p4y - p1y) * t0y;
	if ((t1 >= 0 && t2 <= 0) || (t1 <= 0 && t2 >= 0))
		return true;
	return false;
}

/*
	function to search items in rectangle
*/
#if defined(CALC_CIRCLE) || defined(CALC_POINT)
indexer* search_rect2(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, bool intersection, /*out*/indexer *count_items, ret_callback2_circle callback2, void *data)
#else
indexer* search_rect2(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, bool intersection, /*out*/indexer *count_items)
#endif
{
#ifdef MINIMAL_DEBUG2
	unsigned cc1 = 0, cc2 = 0, cc3 = 0, cc4 = 0;
#endif
	// memory for result
	size_t mem_size = 256;
	size_t count_mem = 1;
	alignas(16) indexer* idxs = (indexer*)aligned_alloc(16, sizeof(indexer) * mem_size * count_mem);
	indexer idx = 0;

	struct node *stack_node[64];
	int stack_pos = 0;
	indexer stack_idx[64];

	coord tx = 0.0, ty = 0.0, tx1 = 0.0, ty1 = 0.0;
	indexer tn = 0, tn1 = 0;

	indexer i = 0;

	coord dist = 0.0;

#ifdef CALC_POINT
	coord tmp_dist = FLT_MAX;
	indexer tmp_idx = -1;
#endif // CALC_POINT
	/*	struct t_result tres;
	coord q1 = nd->x2 - nd->x1;
	coord q2 = nd->y2 - nd->y1;
	tres.dist = sqrt(q1 * q1 + q2 * q2);
	tres.idx = (indexer)-1;
	*/
	while (i < nd->count_child_nodes) {
		// node in bounrary or bounrary in node
		if (nd->x1 <= x_max && nd->x2 >= x_min && nd->y1 <= y_max && nd->y2 >= y_min) {
			// check node fully in the boundary
#if defined(CALC_CIRCLE) || defined(CALC_POINT)
			if (!callback2 && nd->x1 >= x_min && nd->y1 >= y_min && nd->x2 <= x_max && nd->y2 <= y_max) {
#else
			if (nd->x1 >= x_min && nd->y1 >= y_min && nd->x2 <= x_max && nd->y2 <= y_max) {
#endif // CALC_CIRCLE
				// node is fully in the boundary
				add_nodes(nd, &mem_size, &count_mem, &idx, &idxs);
#ifdef MINIMAL_DEBUG2
				printf("Increment node on %u\n", idx - cc3);
				cc3 = idx;
#endif // MINIMAL_DEBUG2
			}
			else {
				// node not fully in the boundaty
				if (nd->is_last_node) {
					for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
						struct branch *br = (struct branch*)(nd->child_node)[j];
						// checking like node
						if (br->x_min <= x_max && br->x_max >= x_min && br->y_min <= y_max && br->y_max >= y_min) {
							// check branch fully in the boundary
#if defined(CALC_CIRCLE) || defined(CALC_POINT)
							if (!callback2 && br->x_min >= x_min && br->y_min >= y_min && br->x_max <= x_max && br->y_max <= y_max) {
#else
							if (br->x_min >= x_min && br->y_min >= y_min && br->x_max <= x_max && br->y_max <= y_max) {
#endif // CALC_CIRCLE
								// branch is fully in the boundary
								add_branch(br, &mem_size, &count_mem, &idx, &idxs);
#ifdef MINIMAL_DEBUG2
								printf("Increment branch on %u\n", idx - cc3);
								cc3 = idx;
#endif // MINIMAL_DEBUG2
							}
							else {
								for (indexer i1 = 0; i1 < br->count_shapes; ++i1) {
									// checking like node
									if (br->xsh_min[i1] <= x_max && br->xsh_max[i1] >= x_min && br->ysh_min[i1] <= y_max && br->ysh_max[i1] >= y_min) {
										// check shape fully in boundary or intersection with boundary
										bool fl1 = false;
										indexer end = i1 + 1 >= br->count_shapes ? br->count_leafs : br->offset[i1 + 1];
										for (indexer k = br->offset[i1]; k < end; ++k) {
											if (br->leaf_x[k] >= x_min && br->leaf_x[k] <= x_max && br->leaf_y[k] >= y_min && br->leaf_y[k] <= y_max) {
												fl1 = true;
												break;
											}
										}
										if (!fl1 && intersection) {
											// last check: intersection

											// side 1/2
											for (indexer k = br->offset[i1]; k < end; ++k) {
												if (k != end - 1) {
													if (check_intersection(br->leaf_x[k], br->leaf_y[k], br->leaf_x[k + 1], br->leaf_y[k + 1], x_min, y_min, x_max, y_min)) {
														fl1 = true;
														break;
													}
												}
												else {
													if (check_intersection(br->leaf_x[k], br->leaf_y[k], br->leaf_x[br->offset[i1]], br->leaf_y[br->offset[i1]], x_min, y_min, x_max, y_min)) {
														fl1 = true;
														break;
													}
												}
											}

											// side 2/3
											if (!fl1) {
												for (indexer k = br->offset[i1]; k < end; ++k) {
													if (k != end - 1) {
														if (check_intersection(br->leaf_x[k], br->leaf_y[k], br->leaf_x[k + 1], br->leaf_y[k + 1], x_max, y_min, x_max, y_max)) {
															fl1 = true;
															break;
														}
													}
													else {
														if (check_intersection(br->leaf_x[k], br->leaf_y[k], br->leaf_x[br->offset[i1]], br->leaf_y[br->offset[i1]], x_max, y_min, x_max, y_max)) {
															fl1 = true;
															break;
														}
													}
												}
											}

											// side 3/4
											if (!fl1) {
												for (indexer k = br->offset[i1]; k < end; ++k) {
													if (k != end - 1) {
														if (check_intersection(br->leaf_x[k], br->leaf_y[k], br->leaf_x[k + 1], br->leaf_y[k + 1], x_max, y_max, x_min, y_max)) {
															fl1 = true;
															break;
														}
													}
													else {
														if (check_intersection(br->leaf_x[k], br->leaf_y[k], br->leaf_x[br->offset[i1]], br->leaf_y[br->offset[i1]], x_max, y_max, x_min, y_max)) {
															fl1 = true;
															break;
														}
													}
												}
											}

											// side 4/1
											if (!fl1) {
												for (indexer k = br->offset[i1]; k < end; ++k) {
													if (k != end - 1) {
														if (check_intersection(br->leaf_x[k], br->leaf_y[k], br->leaf_x[k + 1], br->leaf_y[k + 1], x_min, y_max, x_min, y_min)) {
															fl1 = true;
															break;
														}
													}
													else {
														if (check_intersection(br->leaf_x[k], br->leaf_y[k], br->leaf_x[br->offset[i1]], br->leaf_y[br->offset[i1]], x_min, y_max, x_min, y_min)) {
															fl1 = true;
															break;
														}
													}
												}
											}
										}

#ifdef MINIMAL_DEBUG2
										if (fl1)
											cc4++;
#endif

#ifdef CALC_CIRCLE
										// check callback function to store in return collection
										if (fl1 && callback2 && data) {
											struct calc_data cc;
											cc.br = br;
											cc.idx = i1;
											struct point *center = (struct point*)data;
											cc.center = struct point(*center);
											cc.radius = (x_max - x_min) / 2.0;

											fl1 = callback2(&cc);
#ifdef MINIMAL_DEBUG2
											cc1++;
											if (fl1)
												cc2++;
#endif // MINIMAL_DEBUG2
										}
#elif defined(CALC_POINT)
										// check callback function to store in return collection
										if (fl1 && callback2 && data) {
											struct calc_data cc;
											cc.br = br;
											cc.idx = i1;
											struct point *center = (struct point*)data;
											cc.center = point(*center);
											//cc.radius = (x_max - x_min) / 2.0;
											cc.dist = FLT_MAX;
											cc.curr_idx = -1;

											fl1 = callback2(&cc);
											if (cc.dist < tmp_dist) {
												tmp_dist = cc.dist;
												tmp_idx = cc.curr_idx;
											}
											idx = 1;
										}

#endif //CALC_CIRCLE

										if (fl1) {
#ifdef MINIMAL_DEBUG2
											cc3++;
#endif // MINIMAL_DEBUG2
											// check memory and store index
											if (idx > mem_size * count_mem) {
												count_mem++;
												idxs = (indexer*)_aligned_realloc(idxs, sizeof(indexer) * mem_size * count_mem, 16);
											}
											// store index of current sergment from shape
											idxs[idx] = br->leaf_number[br->offset[i1]];
											idx++;

											/*indexer to_end;
											if (i1 == br->count_shapes - 1)
												to_end = br->count_leafs;
											else
												to_end = br->offset[i1 + 1];
												*/
										}
									}
								}
							}
						}
					}
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
					}
				}
				else {
					stack_idx[stack_pos] = i;
					stack_node[stack_pos] = nd;
					stack_pos++;
					i = 0;
					nd = (struct node*)nd->child_node[i];
				}
			}
					/*} else if (i < nd->count_child_nodes) {
					i++;*/
		}
		else if (stack_pos > 0) {
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
			}
		}
		else {
			char ch1[128];
#ifndef _WIN
			snprintf(ch1, 128, "node: x_min=%.2f, y_min=%.2f, x_max=%.2f, y_max=%.2f, rect: x_min=%.2f, y_min=%.2f, x_max=%.2f, y_max=%.2f", nd->x1, nd->y1, nd->x2, nd->x2, x_min, y_min, x_max, y_max);
#else
			sprintf_s(ch1, 128, "node: x_min=%.2f, y_min=%.2f, x_max=%.2f, y_max=%.2f, rect: x_min=%.2f, y_min=%.2f, x_max=%.2f, y_max=%.2f", nd->x1, nd->y1, nd->x2, nd->x2, x_min, y_min, x_max, y_max);
#endif
			lprintf(ch1);
			return NULL; //(indexer)-1;
		}
		//continue;
	};

	//free(stack_node);
	//free(stack_idx);
#ifdef MINIMAL_DEBUG
	char ch1[128];
	//sprintf_s(ch1, 128, "x=%.2f, y=%.2f, point: x=%.2f, y=%.2f, x=%.2f, y=%.2f dist=%.2f, n1=%u, n2=%u", x, y, tx, ty, tx1, ty1, tres.dist, tn, tn1);
	sprintf_s(ch1, 128, "c1 = %u, c2 = %u, utilization = %0.2f", temp_counter1, temp_counter2, temp_counter2 * 100.0 / temp_counter1);
	lprintf(ch1);
#endif
	*count_items = idx;
	idxs = (indexer*)_aligned_realloc(idxs, sizeof(indexer) * idx, 16);
#ifdef MINIMAL_DEBUG2
	printf("DEBUG: count1 = %u, count2 = %u, count3 = %u, count4 = %u\n", cc1, cc2, cc3, cc4);
#endif

#ifdef CALC_POINT
	idxs[0] = tmp_idx;
	if (data) {
		struct point *center = (struct point*)data;
		center->x = tmp_dist;
	}
#endif
	return idxs; //(indexer)-1;
}

/*
	function to calculate item in circle
*/
/* ret_callback2_circle */ bool scircle2(void *p)
{
	struct calc_data *cc = (struct calc_data*)p;
	indexer end = cc->idx + 1 >= cc->br->count_shapes ? cc->br->count_leafs : cc->br->offset[cc->idx + 1];
	indexer k = cc->br->offset[cc->idx];
	struct point pc(cc->center);
	if (k == end) {
#ifndef CALC_POINT
		if (sqrt((cc->br->leaf_x[k] - pc.x) * (cc->br->leaf_x[k] - pc.x) + (cc->br->leaf_y[k] - pc.y) + (cc->br->leaf_y[k] - pc.y)) <= cc->radius) {
		//if ((cc->br->leaf_x[k] - pc.x) * (cc->br->leaf_x[k] - pc.x) + (cc->br->leaf_y[k] - pc.y) + (cc->br->leaf_y[k] - pc.y) <= cc->radius * cc->radius) {
			return true;
		}
		else {
			return false;
		}
#else
		cc->dist = sqrt((cc->br->leaf_x[k] - pc.x) * (cc->br->leaf_x[k] - pc.x) + (cc->br->leaf_y[k] - pc.y) + (cc->br->leaf_y[k] - pc.y));
		cc->curr_idx = cc->br->leaf_number[cc->br->offset[cc->idx]];
		return false; // because point
#endif
	}

	struct point p0, p1;
	coord t1;
	for ( ; k < end; ++k) {
		p0 = point(cc->br->leaf_x[k], cc->br->leaf_y[k]);
		if (k != end - 1) 
			p1 = point(cc->br->leaf_x[k + 1], cc->br->leaf_y[k + 1]);
		else
			p1 = point(cc->br->leaf_x[cc->br->offset[cc->idx]], cc->br->leaf_y[cc->br->offset[cc->idx]]);
		t1 = distance2(&pc, &p0, &p1);

#ifndef CALC_POINT
		if (t1 <= cc->radius)
			return true;
#else
		if (t1 < cc->dist) {
			cc->dist = t1;
			cc->curr_idx = cc->br->leaf_number[cc->br->offset[cc->idx]];
		}
#endif
	}

	return false;
}

#ifdef CALC_CIRCLE
/*
	search number of items from circle
*/
indexer* search_circle2(/*in*/struct node *nd, /*in*/coord x, /*in*/coord y, /*in*/coord radius, /*in*/bool intersection, /*out*/indexer *count_items)
{
	//indexer count_tmp, count = 0;
	alignas(16) indexer* idxs_tmp = search_rect2(nd, x - radius, y - radius, x + radius, y + radius, intersection, count_items, scircle2, &(struct point(x, y)));
	
	/*__declspec(align(16)) indexer* idxs = (indexer*)_aligned_malloc(sizeof(indexer) * count_tmp, 16);

	struct branch* br = NULL;
	struct node* nd1 = nd;
	while (!nd1->is_last_node) {
		nd1 = (struct node*)(nd1->child_node[0]);
	}
	br = (struct branch*)(nd1->child_node[0]);

	for (indexer i = 0; i < count_tmp; ++i) {
		indexer start = br->offset[idxs_tmp[i]];
		//if ()
	}
	*/
	return idxs_tmp;
}
#elif defined(CALC_POINT)
indexer* search_nearest_item2(/*in*/struct node *nd, /*in*/coord x, /*in*/coord y, /*in*/coord radius, bool intersection, /*out*/indexer *count_items, /*out*/coord *dist)
{
	//indexer count_tmp, count = 0;
	struct point tmp(x, y);
	alignas(16) indexer* idxs_tmp = search_rect2(nd, x - radius, y - radius, x + radius, y + radius, intersection, count_items, scircle2, &tmp);
	*dist = tmp.x;
	return idxs_tmp;
}
#endif // CALC_CIrcle
