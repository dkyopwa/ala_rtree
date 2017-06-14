#include <stdio.h>
#include <malloc.h>
#include <stdbool.h>
#include <math.h>
#include "first.h"

//#define MINIMAL_DEBUG
//#define FIND_V1
//#define FIND_V2
// FIND_V2 on the test_data not found hi performance

struct t_result {
	coord dist;
	indexer idx;
};

struct point{
	coord x;
	coord y;
};

// ------------------------------- SERCHIG --------------------------------------
bool search_callback(indexer idx, void *arg)
{
	return true;
}

/// search rect in selected rect
indexer search_rect(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, ret_callback callback, void *arg)
{
	const size_t mem_size = 128;
	unsigned mem_offset = 1;
	struct node **stack_node = (struct node**)malloc(sizeof(struct node*) * mem_size * mem_offset);
	unsigned stack_pos = 0;
	indexer *stack_idx = (indexer*)malloc(sizeof(indexer) * mem_size * mem_offset);

	//struct node *nd = m_nodes;

	indexer i = 0;
	bool flag_del_branches = false;
	while (i < nd->count_child_nodes) {
		/*if (i >= nd->count_child_nodes) {
		break;
		}*/
		stack_idx[stack_pos] = i;
		stack_node[stack_pos] = nd;
		stack_pos++;
		nd = (struct node*)nd->child_node[i];
		if (nd->is_last_node) {
			for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
				struct branch *br = (struct branch*)(nd->child_node)[i];
				// find

			}
			flag_del_branches = true;

			// return from stack
			stack_pos--;
			nd = stack_node[stack_pos];
			i = stack_idx[stack_pos] + 1;
		}
		else {
			if (flag_del_branches) {
				//free(nd->child_node);
				flag_del_branches = false;
			}
		}
	} while (true);

	return 0;
}

/// calculate distance from point to line
coord distance(struct point *p, struct point *line_p0, struct point *line_p1)
{
	coord vx, vy, wx, wy, c1, c2, b, pbx, pby;

	vx = line_p1->x - line_p0->x;
	vy = line_p1->y - line_p0->y;
	wx = p->x - line_p0->x;
	wy = p->y - line_p0->y;
	c1 = vx * wx + vy * wy;

	if (c1 <= 0) {
		//return sqrt(pow(fabs(p->x - line_p0->x), 2) + pow(fabs(p->y - line_p0->y), 2));
		coord t1 = p->x - line_p0->x;
		coord t2 = p->y - line_p0->y;
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

/// search the nearest item for point
indexer search_point(struct node *nd, coord x, coord y, coord radius)
{
#ifdef MINIMAL_DEBUG
	unsigned temp_counter1 = 0, temp_counter2 = 0;
#endif
	//const size_t mem_size = 64;
	//unsigned mem_offset = 1;
	//struct node **stack_node = (struct node**)malloc(sizeof(struct node*) * mem_size * mem_offset);
	struct node *stack_node[64];
	int stack_pos = 0;
	//indexer *stack_idx = (indexer*)malloc(sizeof(indexer) * mem_size * mem_offset);
	indexer stack_idx[64];

	coord tx = 0.0, ty = 0.0, tx1 = 0.0, ty1 = 0.0;
	indexer tn = 0, tn1 = 0;

	indexer i = 0;
	//bool flag_del_branches = false;

	coord dist = 0.0;
	struct t_result tres;
	//tres.dist = sqrt(pow(nd->x2 - nd->x1, 2) + pow(nd->y2 - nd->y1, 2));
	coord q1 = nd->x2 - nd->x1;
	coord q2 = nd->y2 - nd->y1;
	tres.dist = sqrt(q1 * q1 + q2 * q2);
	tres.idx = (indexer)-1;

	while (i < nd->count_child_nodes) {
#ifdef FIND_V2
		if ((nd->x1 <= x || nd->x1 <= x + radius) && (nd->x2 >= x || nd->x2 >= x - radius) && (nd->y1 <= y || nd->y1 <= y + radius) && (nd->y2 >= y || nd->y2 >= y - radius)) {
#else
		if (nd->x1 <= x + radius && nd->x2 >= x - radius && nd->y1 <= y + radius && nd->y2 >= y - radius) {
#endif // FIND_V2
			if (nd->is_last_node) {
				for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
					struct branch *br = (struct branch*)(nd->child_node)[j];
#ifdef FIND_V2
					if ((br->x_min <= x || br->x_min <= x + radius) && (br->x_max >= x || br->x_max >= x - radius) && (br->y_min <= y || br->y_min <= y + radius) && (br->y_max >= y || br->y_max >= y - radius)) {
#else
					if (br->x_min <= x + radius && br->x_max >= x - radius && br->y_min <= y + radius && br->y_max >= y - radius) {
#endif // FIND_V2
#ifndef FIND_V1
						for (indexer i1 = 0; i1 < br->count_shapes; ++i1) {
#ifdef FIND_V2
							if ((br->xsh_min[i1] <= x || br->xsh_min[i1] <= x + radius) && (br->xsh_max[i1] >= x || br->xsh_max[i1] >= x - radius) && (br->ysh_min[i1] <= y || br->ysh_min[i1] <= y + radius) && (br->ysh_max[i1] >= y || br->ysh_max[i1] >= y - radius)) {
#else // FIND_V2
							if (br->xsh_min[i1] <= x + radius && br->xsh_max[i1] >= x - radius && br->ysh_min[i1] <= y + radius && br->ysh_max[i1] >= y - radius) {
#endif // FIND_V2
								indexer to_end;
								if (i1 == br->count_shapes - 1)
									to_end = br->count_leafs;
								else
									to_end = br->offset[i1 + 1];
								for (indexer k = br->offset[i1]; k < to_end; ++k) {
#else // FIND_V1
						indexer to_end = br->count_leafs - 1;
						if (!br->merge_next_leaf[br->count_leafs - 2])
							to_end = br->count_leafs;
						for (indexer k = 0; k < to_end; ++k) {
#endif // FIND_V1
#ifdef MINIMAL_DEBUG
							temp_counter1++;
#endif
							// before calculation check boundary radius and length of segment
							coord t1 = 0.0;
							if (br->merge_next_leaf[k] && br->length[k] > radius) {
							}
							else {
								coord c1 = br->leafs[k].x - x;
								coord c2 = br->leafs[k].y - y;
								////coord t1 = sqrt(c1 * c1 + c2 * c2);
								t1 = c1 * c1 + c2 * c2;
								if (t1 > radius * radius)
									continue;
							}
#ifdef MINIMAL_DEBUG
							temp_counter2++;
#endif
							//if (((br->leafs[k].x <= x && br->leafs[k].x > x - radius) || (br->leafs[k].x >= x && br->leafs[k].x < x + radius)) && ((br->leafs[k].y <= y && br->leafs[k].y > y - radius) || (br->leafs[k].y >= y && br->leafs[k].y < y + radius)))
							//	continue;

							// find minimal distance
							if (br->merge_next_leaf[k]) {
								// line
								struct point p, line_p0, line_p1;
								p.x = x;
								p.y = y;
								line_p0.x = br->leafs[k].x;
								line_p0.y = br->leafs[k].y;
								line_p1.x = br->leafs[k + 1].x;
								line_p1.y = br->leafs[k + 1].y;
								dist = distance(&p, &line_p0, &line_p1);
							}
							else {
								// point
								//dist = sqrt(pow(fabs(br->leafs[k].x - x), 2) + pow(fabs(br->leafs[k].y - y), 2));
								/* coord c1 = br->leafs[k].x - x;
								coord c2 = br->leafs[k].y - y;
								dist = sqrt(c1 * c1 + c2 * c2);
								*/
								dist = sqrt(t1);
							}
							if (dist < tres.dist) {
								tres.dist = dist;
								tres.idx = br->leafs[k].number;
								tx = br->leafs[k].x;
								ty = br->leafs[k].y;
								tx1 = br->leafs[k + 1].x;
								ty1 = br->leafs[k + 1].y;
								tn = br->leafs[k].number;
								tn1 = br->leafs[k + 1].number;
							}
							}
#ifndef FIND_V1
							}
						}
#else
					//}
#endif
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
			} else {
				stack_idx[stack_pos] = i;
				stack_node[stack_pos] = nd;
				stack_pos++;
				i = 0;
				nd = (struct node*)nd->child_node[i];
			}
		/*} else if (i < nd->count_child_nodes) {
			i++;*/
		} else if (stack_pos > 0) {
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
		} else {
			char ch1[128];
#ifndef _WIN
			snprintf(ch1, 128, "x=%.2f, y=%.2f, point: x=%.2f, y=%.2f, x=%.2f, y=%.2f dist=%.2f, n1=%u, n2=%u", x, y, tx, ty, tx1, ty1, tres.dist, tn, tn1);
#else
			sprintf_s(ch1, 128, "x=%.2f, y=%.2f, point: x=%.2f, y=%.2f, x=%.2f, y=%.2f dist=%.2f, n1=%u, n2=%u", x, y, tx, ty, tx1, ty1, tres.dist, tn, tn1);
#endif
			lprintf(ch1);
			return tres.idx; //(indexer)-1;
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
	return tres.idx; //(indexer)-1;
}
