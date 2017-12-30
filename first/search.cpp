#include <stdio.h>
//#include <malloc.h>
#include <stdbool.h>
#include <math.h>
#include "unimem.h"
#include "first.h"

#include <xmmintrin.h>
#ifdef _WIN
#include <intrin.h>
#endif

//#define MINIMAL_DEBUG
//#define FIND_V1
//#define FIND_V2
// FIND_V2 on the test_data not found hi performance

extern "C" coord distance_sse_v4(__m128 *vec1, __m128 *vec2);
extern "C" coord distance_sse_v5(__m128 *vec1, __m128 *vec2, __m128 *tmp);
extern "C" coord distance_sse_v6(__m128 *vec1, __m128 *vec2);

struct t_result {
	coord dist;
	indexer idx;
};

struct point {
	coord x;
	coord y;
};

struct points_distanse_sse {
	__m128 point_x;
	__m128 point_y;
	__m128 line_p0_x;
	__m128 line_p0_y;
	__m128 line_p1_x;
	__m128 line_p1_y;
};

/*alignas(16) struct vector {
	float f1;
	float f2;
	float f3;
	float f4;
};*/
// ------------------------------- SERCHIG --------------------------------------
/* bool search_callback(indexer idx, void *arg)
{
	return true;
}
*/

/// search rect in selected rect
indexer search_rect(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max/*, ret_callback callback, void *arg*/)
{
	const size_t mem_size = 128;
	unsigned mem_offset = 1;
	alignas(16) struct node **stack_node = (struct node**)aligned_alloc(16, sizeof(struct node*) * mem_size * mem_offset);
	unsigned stack_pos = 0;
	alignas(16) indexer *stack_idx = (indexer*)aligned_alloc(16, sizeof(indexer) * mem_size * mem_offset);

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
			/*for (unsigned i = 0; i < nd->count_child_nodes; ++i) {
				struct branch *br = (struct branch*)(nd->child_node)[i];
				// find

			}*/
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
		//coord t1 = p->x - line_p0->x;
		//coord t2 = p->y - line_p0->y;
		coord t1 = wx;
		coord t2 = wy;
		return (coord)sqrt(t1 * t1 + t2 * t2);
	}

	c2 = vx * vx + vy * vy;
	if (c2 <= c1) {
		//return sqrt(pow(fabs(p->x - line_p1->x), 2) + pow(fabs(p->y - line_p1->y), 2));
		coord t1 = p->x - line_p1->x;
		coord t2 = p->y - line_p1->y;
		return (coord)sqrt(t1 * t1 + t2 * t2);
	}

	b = c1 / c2;
	pbx = line_p0->x + b * vx;
	pby = line_p0->y + b * vy;

	//return sqrt(pow(fabs(p->x - pbx), 2) + pow(fabs(p->y - pby), 2));
	coord t1 = p->x - pbx;
	coord t2 = p->y - pby;
	return (coord)sqrt(t1 * t1 + t2 * t2);
}

/// calculate distance from point to line
//coord distance_sse3(struct point *p, struct point *line_p0, struct point *line_p1)
coord distance_sse_v3(__m128 *vec1, __m128 *vec2) // vec1 = p1.x, p1.y, p.x, p.x; vec2 = p0.x, p0.y, p0.x, p0.y
{
	__m128 tmp[3];
#ifdef _WIN
	tmp[0] = {1.0, 2.0, 3.0, 4.0};
	tmp[1] = { 1.0, 2.0, 3.0, 4.0 };
	tmp[2] = { 1.0, 2.0, 3.0, 4.0 };
	coord t3 = distance_sse_v5(vec1, vec2, tmp);
#else
	float *f1; f1 = (float*)&(tmp[0]);
	f1[0] = 1.0; f1[1] = 2.0; f1[2] = 3.0; f1[3] = 4.0;
	float *f2; f2 = (float*)&(tmp[1]);
	f2[0] = 1.0; f2[1] = 2.0; f2[2] = 3.0; f2[3] = 4.0 ;
	float *f3; f3 = (float*)&(tmp[2]);
	f3[0] = 1.0; f3[1] = 2.0; f3[2] = 3.0; f3[3] = 4.0;
	coord t3 = 0.0;
#endif // _WIN

	//coord t1 = distance_sse_v4(vec1, vec2);
	register __m128 res2 = _mm_movehl_ps(*vec1, *vec1); // px, py, px, py
	/*coord q1 = ((float*)&res2)[0];
	coord q2 = ((float*)&res2)[1];
	coord q3 = q1 + q2;
	*/
	register __m128 res1 = _mm_sub_ps(*vec1, *vec2); // where 1 - vx, 2 - vy, 3 - wx(t1(c1)), 4 - wy(t2(c1))
	__m128 res3 = _mm_movelh_ps(res1, res1); // vx, vy, vx, vy
	register __m128 res4 = _mm_sub_ps(res2, *vec1); // t1(c2), t2(c2), 0(unk), 0(unk)
	register __m128 res5 = _mm_mul_ps(res3, res1); // vx*vx, vy*vy, wx*vx, wy*vy
	__m128 res8 = _mm_shuffle_ps(res1, res4, 78); // for sqrt: t1(c1), t2(c1), t1(c2), t2(c2)
	__m128 res6 = _mm_shuffle_ps(res5, res5, 245); // vy*vy, vy*vy, wy*wy, wy*wy
	__m128 res9 = _mm_mul_ps(res8, res8); // for sqrt: t1(c1)^2, t2(c1)^2, t1(c2)^2, t2(c2)^2
	__m128 res7 = _mm_add_ps(res5, res6); // c1, unk, c2, unk
	__m128 res10 = _mm_shuffle_ps(res9, res9, 245); // t2(c1)^2, t2(c1)^2, t2(c2)^2, t2(c2)^2
	__m128 res11 = _mm_add_ps(res9, res10); // t1 * t1 + t2 * t2 (c1), unk, t1 * t1 + t2 * t2 (c2), unk
	__m128 res12 = _mm_sqrt_ps(res11); // sqrt(c1), unk, sqrt(c2), unk

	if (*(float*)&res7 <= 0) {
		return *(float*)&res12;
	}

	if (((float*)&res7)[2] <= *(float*)&res7) {
		return ((float*)&res12)[2];
	}

	//coord b = c1 / c2;
	coord b = *(float*)&res7 / ((float*)&res7)[2];
	__m128 bb = _mm_set1_ps(b);
	//pbx = line_p0->x + b * vx;
	//pby = line_p0->y + b * vy;
	res9 = _mm_mul_ps(bb, res1); // vx * b, vy * b, unk, unk
	res10 = _mm_add_ps(*vec2, res9); // pbx, pby, unk, unk

	//coord t1 = p->x - pbx;
	//coord t2 = p->y - pby;
	res11 = _mm_sub_ps(res2, res10); // t1, t2, unk, unk
	res12 = _mm_mul_ps(res11, res11); // t1^2, t2^2, unk, unk
	/* __m128 res13 = _mm_shuffle_ps(res12, res12, 225);
	__m128 res14 = _mm_add_ps(res13, res12);
	__m128 res15 = _mm_sqrt_ss(res14);
	//return sqrt(t1 * t1 + t2 * t2);
	return *(float*)&res15; */
	coord t2 = (coord)sqrt(*((float*)&res12) + ((float*)&res12)[1]);
	if (t3 != t2)
		return t2;
	return t2;
}

/// calculate distance from point to line
//coord distance_sse2(struct point *p, struct point *line_p0, struct point *line_p1)
coord distance_sse_v2(__m128 *vec1, __m128 *vec2) // vec1 = p1.x, p1.y, p.x, p.x; vec2 = p0.x, p0.y, p0.x, p0.y
{
	//coord vx, vy, wx, wy, c1, c2, b, pbx, pby;

	__m128 res1 = _mm_sub_ps(*vec1, *vec2); // where 1 - vx, 2 - vy, 3 - wx, 4 - wy
	__m128 tt1 = _mm_movehl_ps(res1, res1); // wx, wy, wx, wy
	__m128 res2 = _mm_mul_ps(res1, tt1); // c1.1, c1.2, ret_c1 t1*t1, ret_c1 t2*t2

/*	vx = line_p1->x - line_p0->x;
	vy = line_p1->y - line_p0->y;
	wx = p->x - line_p0->x;
	wy = p->y - line_p0->y;
	c1 = vx * wx + vy * wy;
*/
	__m128 res3 = _mm_shuffle_ps(res2, res2, 177); // c1.2, c1.1, ret_c1 t2*t2, ret_c1 t1*t1
	__m128 res4 = _mm_add_ps(res2, res3); // c1, c1, for sqrt, for sqrt

	//if (c1 <= 0) {
	if ((*((float*)&res4)) <= 0) {
		//coord t1 = p->x - line_p0->x;
		//coord t2 = p->y - line_p0->y;
/*		coord t1 = wx;
		coord t2 = wy;
		return sqrt(t1 * t1 + t2 * t2); */
		coord t1 = (coord)sqrt(((float*)&res4)[3]);
	/*	coord t2 = distance_sse_v6(vec1, vec2);
		if (t1 != t2)
			return t1;*/
		return t1;
	}

	//c2 = vx * vx + vy * vy;
	__m128 res5 = _mm_mul_ps(res1, res1); // vx^2, vy^2, unk, unk
	__m128 res51 = _mm_shuffle_ps(res5, res5, 225); // vy^2, vx^2, unk, unk
	__m128 res52 = _mm_add_ss(res5, res51); // c2, unk, unk, unk
	//__m128 res6 = _mm_setzero_ps();
	__m128 res6 = _mm_shuffle_ps(*vec1, *vec1, 78); // p.x, p.y, p1.x, p1.y
	//if (c2 <= c1) {
	if (*((float*)&res52) <= *((float*)&res4)) {
		//coord t1 = p->x - line_p1->x;
		//coord t2 = p->y - line_p1->y;
		//res6 = _mm_shuffle_ps(*vec1, *vec1, 78); // p.x, p.y, p1.x, p1.y
		__m128 res7 = _mm_sub_ps(res6, res1); // t1, t2, unk, unk
		__m128 res8 = _mm_mul_ps(res7, res7);
		/*__m128 res81 = _mm_shuffle_ps(res7, res7, 225);
		__m128 res82 = _mm_add_ss(res8, res81); */
		//return sqrt(t1 * t1 + t2 * t2);
		coord t1 = (coord)sqrt(*((float*)&res8) + ((float*)&res8)[1]);
	/*	coord t2 = distance_sse_v6(vec1, vec2);
		if (t1 != t2)
			return t1;*/
		return t1;
		/*__m128 res83 = _mm_sqrt_ss(res82);
		return *(float*)&res83; */
	}

	//coord b = c1 / c2;
	coord b = *((float*)&res4) / *((float*)&res52);
	__m128 bb = _mm_set1_ps(b);
	//pbx = line_p0->x + b * vx;
	//pby = line_p0->y + b * vy;
	__m128 res9 = _mm_mul_ps(bb, res1); // vx * b, vy * b, unk, unk
	__m128 res10 = _mm_add_ps(*vec2, res9); // pbx, pby, unk, unk

	//coord t1 = p->x - pbx;
	//coord t2 = p->y - pby;
	__m128 res11 = _mm_sub_ps(res6, res10); // t1, t2, unk, unk
	__m128 res12 = _mm_mul_ps(res11, res11); // t1^2, t2^2, unk, unk
	/* __m128 res13 = _mm_shuffle_ps(res12, res12, 225);
	__m128 res14 = _mm_add_ps(res13, res12);
	__m128 res15 = _mm_sqrt_ss(res14);
	//return sqrt(t1 * t1 + t2 * t2);
	return *(float*)&res15; */
	coord t1 = (coord)sqrt(*((float*)&res12) + ((float*)&res12)[1]);
/*	coord t2 = distance_sse_v6(vec1, vec2);
	if (t1 != t2)
		return t1;*/

	return t1;
}

/// calculate distance from point to line
__m128 distance_sse(points_distanse_sse *data)
{
	__m128 vx, vy, wx, wy, c1, c2, b, pbx, pby;
	__m128 t1, t2;
	__m128 result;

	//vx = data->line_p1_x - data->line_p0_x;
	vx = _mm_sub_ps(data->line_p1_x, data->line_p0_x);
	//vy = data->line_p1_y - data->line_p0_y;
	vy = _mm_sub_ps(data->line_p1_y, data->line_p0_y);
	//wx = data->point_x - data->line_p0_x;
	wx = _mm_sub_ps(data->point_x, data->line_p0_x);
	//wy = data->point_y - data->line_p0_y;
	wy = _mm_sub_ps(data->point_y, data->line_p0_y);
	//c1 = vx * wx + vy * wy;
	c1 = _mm_add_ps(_mm_mul_ps(vx, wx), _mm_mul_ps(vy, wy));

	/*if (c1 <= 0) {
		//return sqrt(pow(fabs(p->x - line_p0->x), 2) + pow(fabs(p->y - line_p0->y), 2));
		coord t1 = p->x - line_p0->x;
		coord t2 = p->y - line_p0->y;
		return sqrt(t1 * t1 + t2 * t2);
	}*/
	__m128 tmp_c1 = _mm_setzero_ps();
	__m128 res_c1 = _mm_cmple_ps(c1, tmp_c1);
	///t1 = _mm_sub_ps(data->point_x, data->line_p0_x);
	///t2 = _mm_sub_ps(data->point_y, data->line_p0_y);
	///result = _mm_rsqrt_ps(_mm_add_ps(_mm_mul_ps(t1, t1), _mm_mul_ps(t2, t2)));
	result = _mm_sqrt_ps(_mm_add_ps(_mm_mul_ps(wx, wx), _mm_mul_ps(wy, wy)));
	unsigned *temp1 = reinterpret_cast<unsigned*>(&res_c1);
	if (temp1[0] && temp1[1] && temp1[2] && temp1[3]) {
		return result;
	}
	
/*	if (*temp1 && *(temp1 + 1) && *(temp1 + 2) && *(temp1 + 3)) {
		return result;
	} */

	/*c2 = vx * vx + vy * vy;
	if (c2 <= c1) {
		//return sqrt(pow(fabs(p->x - line_p1->x), 2) + pow(fabs(p->y - line_p1->y), 2));
		coord t1 = p->x - line_p1->x;
		coord t2 = p->y - line_p1->y;
		return sqrt(t1 * t1 + t2 * t2);
	}
	*/
	c2 = _mm_add_ps(_mm_mul_ps(vx, vx), _mm_mul_ps(vy, vy));
	t1 = _mm_sub_ps(data->point_x, data->line_p1_x);
	t2 = _mm_sub_ps(data->point_y, data->line_p1_y);
	__m128 result2 = _mm_sqrt_ps(_mm_add_ps(_mm_mul_ps(t1, t1), _mm_mul_ps(t2, t2)));
	__m128 res_c2 = _mm_cmple_ps(c2, c1);
	unsigned *temp2 = reinterpret_cast<unsigned*>(&res_c2);

	//b = c1 / c2;
	b = _mm_div_ps(c1, c2);
	//pbx = line_p0->x + b * vx;
	pbx = _mm_add_ps(data->line_p0_x, _mm_mul_ps(b, vx));
	//pby = line_p0->y + b * vy;
	pby = _mm_add_ps(data->line_p0_y, _mm_mul_ps(b, vy));

	//return sqrt(pow(fabs(p->x - pbx), 2) + pow(fabs(p->y - pby), 2));
	//coord t1 = p->x - pbx;
	t1 = _mm_sub_ps(data->point_x, pbx);
	//coord t2 = p->y - pby;
	t2 = _mm_sub_ps(data->point_y, pby);
	//return sqrt(t1 * t1 + t2 * t2);
	__m128 result3 = _mm_sqrt_ps(_mm_add_ps(_mm_mul_ps(t1, t1), _mm_mul_ps(t2, t2)));

/*	unsigned w1[2]; w1[0] = 0, w1[1] = 1;
	if (w1[1])
		if (w1[0])
			w1[0] = w1[1];
		else
			w1[1] = w1[0];
*/
	if (!temp1[0]) {
		if (temp2[0])
			((float*)&result)[0] = ((float*)&result2)[0];
		else
			((float*)&result)[0] = ((float*)&result3)[0];
	}
	if (!temp1[1]) {
		if (temp2[1])
			((float*)&result)[1] = ((float*)&result2)[1];
		else
			((float*)&result)[1] = ((float*)&result3)[1];
	}
	if (!temp1[2]) {
		if (temp2[2])
			((float*)&result)[2] = ((float*)&result2)[2];
		else
			((float*)&result)[2] = ((float*)&result3)[2];
	}
	if (!temp1[3]) {
		if (temp2[3])
			((float*)&result)[3] = ((float*)&result2)[3];
		else
			((float*)&result)[3] = ((float*)&result3)[3];
	}

	return result;
}

/// search the nearest item for point
indexer search_point_sse(struct node *nd, coord x, coord y, coord radius)
{
#ifdef MINIMAL_DEBUG
	unsigned temp_counter1 = 0, temp_counter2 = 0;
#endif

/*	__m128 px = _mm_set1_ps(5.0);
	__m128 py = _mm_set1_ps(5.0);
	__m128 line0_x = { 1.0f, 3.0f, 5.0f, 10.0f };
	__m128 line0_y = { 2.0f, 2.0f, 2.0f, 2.0f };
	__m128 line1_x = { 7.0f, 10.0f, 5.0f, 7.0f };
	__m128 line1_y = { 20.0f, 10.0f, 7.0f, 10.0f };
	struct points_distanse_sse poi;
	poi.point_x = px;
	poi.point_y = py;
	poi.line_p0_x = line0_x;
	poi.line_p0_y = line0_y;
	poi.line_p1_x = line1_x;
	poi.line_p1_y = line1_y;
	unsigned __int64 t1 = __rdtsc();
	__m128 res1 = distance_sse(&poi);
	unsigned __int64 t2 = __rdtsc();

	struct point p; p.x = 5; p.y = 5;
	struct point p1[4];
	p1[0].x = 1; p1[1].x = 3; p1[2].x = 5; p1[3].x = 10;
	p1[0].y = 2; p1[1].y = 2; p1[2].y = 2; p1[3].y = 2;
	struct point p2[4];
	p2[0].x = 7; p2[1].x = 10; p2[2].x = 5; p2[3].x = 7;
	p2[0].y = 20; p2[1].y = 10; p2[2].y = 7; p2[3].y = 10;
	float res2[4];
	__int64 t3 = __rdtsc();
	for (unsigned i = 0; i < 4; ++i) {
		res2[i] = distance(&p, &(p1[i]), &(p2[i]));
	}
	__int64 t4 = __rdtsc();

	float res3[4];
	__int64 t5 = __rdtsc();
	for (unsigned i = 0; i < 4; ++i) {
		__m128 vec1, vec2;
		//((float*)&vec1)[0] = p2[i].x;
		//((float*)&vec1)[1] = p2[i].y;
		//((float*)&vec1)[2] = p.x;
		//((float*)&vec1)[3] = p.y;
		//((float*)&vec2)[0] = p1[i].x;
		//((float*)&vec2)[1] = p1[i].y;
		//((float*)&vec2)[2] = p1[i].x;
		//((float*)&vec2)[3] = p1[i].y;
		
		vec1 = _mm_setr_ps(p2[i].x, p2[i].y, p.x, p.y);
		vec2 = _mm_setr_ps(p1[i].x, p1[i].y, p1[i].x, p1[i].y);
		res3[i] = distance_sse_v2(&vec1, &vec2);
	}
	__int64 t6 = __rdtsc();

	t1 = t2 - t1;
	t3 = t4 - t3;
	t5 = t6 - t5;
	return -1;
*/
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
	tres.dist = (coord)sqrt(q1 * q1 + q2 * q2);
	tres.idx = (indexer)-1;

	while (i < nd->count_child_nodes) {
		if (nd->x1 <= x + radius && nd->x2 >= x - radius && nd->y1 <= y + radius && nd->y2 >= y - radius) {
			if (nd->is_last_node) {
				for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
					struct branch *br = (struct branch*)(nd->child_node)[j];
					if (br->x_min <= x + radius && br->x_max >= x - radius && br->y_min <= y + radius && br->y_max >= y - radius) {
						for (indexer i1 = 0; i1 < br->count_shapes; ++i1) {
							if (br->xsh_min[i1] <= x + radius && br->xsh_max[i1] >= x - radius && br->ysh_min[i1] <= y + radius && br->ysh_max[i1] >= y - radius) {
								indexer to_end;
								if (i1 == br->count_shapes - 1)
									to_end = br->count_leafs;
								else
									to_end = br->offset[i1 + 1];
								for (indexer k = br->offset[i1]; k < to_end; ++k) {
#ifdef MINIMAL_DEBUG
									temp_counter1++;
#endif
									// before calculation check boundary radius and length of segment
									coord t1 = 0.0;
									if (br->merge_next_leaf[k] && br->length[k] > radius) {
									} else {
										// TO DO LEAFS
										coord c1 = br->leaf_x[k] - x;
										coord c2 = br->leaf_y[k] - y;
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
										// TO DO LEAFS
										line_p0.x = br->leaf_x[k];
										line_p0.y = br->leaf_y[k];
										line_p1.x = br->leaf_x[k + 1];
										line_p1.y = br->leaf_y[k + 1];
										dist = distance(&p, &line_p0, &line_p1);
										
										/*struct vector v1, v2;
										v1.f1 = br->leaf_x[k + 1]; v1.f2 = br->leaf_y[k + 1]; v1.f3 = x; v1.f4 = y;
										v2.f1 = br->leaf_x[k]; v2.f2 = br->leaf_y[k]; v2.f3 = br->leaf_x[k]; v2.f4 = br->leaf_y[k];
										__m128 vec1 = _mm_load_ps((const float*)&v1); // _mm_setr_ps(br->leaf_x[k + 1], br->leaf_y[k + 1], x, y);
										__m128 vec2 = _mm_load_ps((const float*)&v2); // _mm_setr_ps(br->leaf_x[k], br->leaf_y[k], br->leaf_x[k], br->leaf_y[k]);
										dist = distance_sse_v6(&vec1, &vec2); */
									}
									else {
										// point
										//dist = sqrt(pow(fabs(br->leafs[k].x - x), 2) + pow(fabs(br->leafs[k].y - y), 2));
										/* coord c1 = br->leafs[k].x - x;
										coord c2 = br->leafs[k].y - y;
										dist = sqrt(c1 * c1 + c2 * c2);
										*/
										dist = (coord)sqrt(t1);
									}
									if (dist < tres.dist) {
										tres.dist = dist;
										// TO DO LEAFS
										tres.idx = br->leaf_number[k];
										tx = br->leaf_x[k];
										ty = br->leaf_y[k];
										tx1 = br->leaf_x[k + 1];
										ty1 = br->leaf_y[k + 1];
										tn = br->leaf_number[k];
										tn1 = br->leaf_number[k + 1];
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

/// search the nearest item for point
indexer search_point(struct node *nd, coord x, coord y, coord radius)
{
#ifdef MINIMAL_DEBUG
	unsigned temp_counter1 = 0, temp_counter2 = 0;
#endif
	return search_point_sse(nd, x, y, radius);
//	return -1;

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
	tres.dist = (coord)sqrt(q1 * q1 + q2 * q2);
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
#ifdef OLD_LEAFS
								coord c1 = br->leafs[k].x - x;
								coord c2 = br->leafs[k].y - y;
#else
								// TO DO LAEFS
								coord c1 = br->leaf_x[k] - x;
								coord c2 = br->leaf_y[k] - y;
#endif // OLD_LEAFS
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
#ifdef OLD_LEAFS
								line_p0.x = br->leafs[k].x;
								line_p0.y = br->leafs[k].y;
								line_p1.x = br->leafs[k + 1].x;
								line_p1.y = br->leafs[k + 1].y;
#else
								// TO DO LEAFS
								line_p0.x = br->leaf_x[k];
								line_p0.y = br->leaf_y[k];
								line_p1.x = br->leaf_x[k + 1];
								line_p1.y = br->leaf_y[k + 1];
#endif // OLD_LEAFS
								dist = distance(&p, &line_p0, &line_p1);
							}
							else {
								// point
								//dist = sqrt(pow(fabs(br->leafs[k].x - x), 2) + pow(fabs(br->leafs[k].y - y), 2));
								/* coord c1 = br->leafs[k].x - x;
								coord c2 = br->leafs[k].y - y;
								dist = sqrt(c1 * c1 + c2 * c2);
								*/
								dist = (coord)sqrt(t1);
							}
							if (dist < tres.dist) {
								tres.dist = dist;
#ifdef OLD_LEAFS
								tres.idx = br->leafs[k].number;
								tx = br->leafs[k].x;
								ty = br->leafs[k].y;
								tx1 = br->leafs[k + 1].x;
								ty1 = br->leafs[k + 1].y;
								tn = br->leafs[k].number;
								tn1 = br->leafs[k + 1].number;
#else
								// TO DO LAEFS
								tres.idx = br->leaf_number[k];
								tx = br->leaf_x[k];
								ty = br->leaf_y[k];
								tx1 = br->leaf_x[k + 1];
								ty1 = br->leaf_y[k + 1];
								tn = br->leaf_number[k];
								tn1 = br->leaf_number[k + 1];
#endif // OLD_LEAFS
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
					}
			else {
				stack_idx[stack_pos] = i;
				stack_node[stack_pos] = nd;
				stack_pos++;
				i = 0;
				nd = (struct node*)nd->child_node[i];
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

/// search items in rectangle (allocated memory have to free after using)
indexer* search_in_rect(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, indexer *count_items)
{
	// memory for result
	size_t mem_size = 256;
	size_t count_mem = 1;
	alignas(16) indexer* idxs = (indexer*)aligned_alloc(16, sizeof(indexer) * mem_size * count_mem);
	indexer idx = 0;

#ifdef MINIMAL_DEBUG
	unsigned temp_counter1 = 0, temp_counter2 = 0;
#endif
	struct node *stack_node[64];
	int stack_pos = 0;
	indexer stack_idx[64];

	//coord tx = 0.0, ty = 0.0, tx1 = 0.0, ty1 = 0.0;
	//indexer tn = 0, tn1 = 0;

	indexer i = 0;

/*	coord dist = 0.0;
	struct t_result tres;
	coord q1 = nd->x2 - nd->x1;
	coord q2 = nd->y2 - nd->y1;
	tres.dist = sqrt(q1 * q1 + q2 * q2);
	tres.idx = (indexer)-1;
*/
	while (i < nd->count_child_nodes) {
		// node in bounrary or bounrary in node
		if (nd->x1 <= x_max && nd->x2 >= x_min && nd->y1 <= y_max && nd->y2 >= y_min) {
			if (nd->is_last_node) {
				for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
					struct branch *br = (struct branch*)(nd->child_node)[j];
					// checking like node
					if (br->x_min <= x_max && br->x_max >= x_min && br->y_min <= y_max && br->y_max >= y_min) {
						for (indexer i1 = 0; i1 < br->count_shapes; ++i1) {
							// checking like node
							if (br->xsh_min[i1] <= x_max && br->xsh_max[i1] >= x_min && br->ysh_min[i1] <= y_max && br->ysh_max[i1] >= y_min) {
								// check memory and store index
								if (idx > mem_size * count_mem) {
									count_mem++;
									idxs = (indexer*)_aligned_realloc(idxs, sizeof(indexer) * mem_size * count_mem, 16);
								}
								// store index of current sergment from shape
#ifdef OLD_LEAFS
								idxs[idx] = br->leafs[br->offset[i1]].number;
#else
								// TO DO LEAFS
								idxs[idx] = br->leaf_number[br->offset[i1]];
#endif // OLD_LEAFS
								idx++;

								/*indexer to_end;
								if (i1 == br->count_shapes - 1)
									to_end = br->count_leafs;
								else
									to_end = br->offset[i1 + 1];*/
#ifdef MINIMAL_DEBUG
								temp_counter1++;
#endif

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
								//printf("stack_pos 1 = %u\n", stack_pos);
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
						//printf("stack_pos 2 = %u\n", stack_pos);
						i = 0;
						nd = (struct node*)nd->child_node[i];
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
						//printf("stack_pos 3 = %u\n", stack_pos);
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

	return idxs; //(indexer)-1;
}

indexer* search_in_circles(/*in*/struct node *nd, /*in*/coord x, /*in*/coord y, /*in*/coord radius, /*out*/indexer *count_items)
{
#ifdef MINIMAL_DEBUG
	unsigned temp_counter1 = 0, temp_counter2 = 0;
#endif
	// memory for result
	size_t mem_size = 256;
	size_t count_mem = 1;
	alignas(16) indexer* idxs = (indexer*)aligned_alloc(16, sizeof(indexer) * mem_size * count_mem);
	indexer idx = 0;

	struct node *stack_node[64];
	int stack_pos = 0;
	indexer stack_idx[64];

	//coord tx = 0.0, ty = 0.0, tx1 = 0.0, ty1 = 0.0;
	coord c1, c2, t1;
	//indexer tn = 0, tn1 = 0;

	indexer i = 0;

	//coord dist = 0.0;
	coord radius2 = radius * radius;
/*	struct t_result tres;
	coord q1 = nd->x2 - nd->x1;
	coord q2 = nd->y2 - nd->y1;
	tres.dist = sqrt(q1 * q1 + q2 * q2);
	tres.idx = (indexer)-1;
*/
	while (i < nd->count_child_nodes) {
		if (nd->x1 <= x + radius && nd->x2 >= x - radius && nd->y1 <= y + radius && nd->y2 >= y - radius) {
			if (nd->is_last_node) {
				for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
					struct branch *br = (struct branch*)(nd->child_node)[j];
					if (br->x_min <= x + radius && br->x_max >= x - radius && br->y_min <= y + radius && br->y_max >= y - radius) {
						for (indexer k = 0; k < br->count_leafs; ++k) {
#ifdef MINIMAL_DEBUG
							temp_counter1++;
#endif
#ifdef OLD_LEAFS
							c1 = br->leafs[k].x - x;
							c2 = br->leafs[k].y - y;
#else
							// TO DO LEAFS
							c1 = br->leaf_x[k] - x;
							c2 = br->leaf_y[k] - y;
#endif //OLD_LEAFS
							////coord t1 = sqrt(c1 * c1 + c2 * c2);
							t1 = c1 * c1 + c2 * c2;
							if (t1 > radius2)
								continue;
							else {
								// check memory and store index
								if (idx > mem_size * count_mem) {
									count_mem++;
									idxs = (indexer*)_aligned_realloc(idxs, sizeof(indexer) * mem_size * count_mem, 16);
								}
								// store index of current sergment from shape
#ifdef OLD_LEAFS
								idxs[idx] = br->leafs[k].number;
#else
								// TO DO LEAFS
								idxs[idx] = br->leaf_number[k];
#endif // OLD_LEAFS

								do {
									++k;
#ifdef OLD_LEAFS
								} while (idxs[idx] == br->leafs[k].number);
#else
									// TO DO LEAFS
								} while (idxs[idx] == br->leaf_number[k]);
#endif // OLD_LRAFS
								idx++;
							}
#ifdef MINIMAL_DEBUG
							temp_counter2++;
#endif
							//if (((br->leafs[k].x <= x && br->leafs[k].x > x - radius) || (br->leafs[k].x >= x && br->leafs[k].x < x + radius)) && ((br->leafs[k].y <= y && br->leafs[k].y > y - radius) || (br->leafs[k].y >= y && br->leafs[k].y < y + radius)))
							//	continue;

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
			//char ch1[128];
#ifndef _WIN
			//snprintf(ch1, 128, "x=%.2f, y=%.2f, point: x=%.2f, y=%.2f, x=%.2f, y=%.2f dist=%.2f, n1=%u, n2=%u", x, y, tx, ty, tx1, ty1, tres.dist, tn, tn1);
#else
			//sprintf_s(ch1, 128, "x=%.2f, y=%.2f, point: x=%.2f, y=%.2f, x=%.2f, y=%.2f dist=%.2f, n1=%u, n2=%u", x, y, tx, ty, tx1, ty1, tres.dist, tn, tn1);
#endif
			//lprintf(ch1);
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

	return idxs;
}
