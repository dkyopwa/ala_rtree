/*
	Project present like RTree
	Created by Vladimir Nedved 2017
	Apache License 2.0
*/
#pragma once

#ifndef SORT_H_HEADERS
#define SORT_H_HEADERS

#include "first.h"
#include "pch.h"

#ifndef _WIN
// only for 8-bits digit
inline void min_max_sort(int &a, int &b) // https://habrahabr.ru/post/204682/
{
    int d = a - b;
    int m = ~(d >> 8);
    b += d&m;
    a -= d&m;
}

inline void min_max_sort(coord &a, coord &b) // сортирует пару чисел
{
    if(a > b)
    {
        coord t = a;
        a = b;
        b = t;
    }
}
#else 
/// only for 8-bits digit
/* inline void min_max_sort(int *a, int *b) // https://habrahabr.ru/post/204682/
{
	int d = *a - *b;
	int m = ~(d >> 8);
	*b += d&m;
	*a -= d&m;
}

inline void min_max_sort(coord *a, coord *b) // сортирует пару чисел
{
	if (*a > *b)
	{
		coord t = *a;
		*a = *b;
		*b = t;
	}
}
*/
#endif // _WIN

#define MAX_STACK 64 // max size stack
template<class T>
void q_sort_i(T a[], long size) {

	long i, j; // pointers to separated parts
	long lb, ub; // boundary of sorted fragment

	long lbstack[MAX_STACK], ubstack[MAX_STACK]; // стек запросов
											   // каждый запрос задается парой значений,
											   // а именно: левой(lbstack) и правой(ubstack)
											   // границами промежутка
	long stackpos = 1; // текущая позиция стека
	long ppos; // середина массива
	T pivot; // опорный элемент
	T temp;

	lbstack[1] = 0;
	ubstack[1] = size - 1;

	do {
		// Взять границы lb и ub текущего массива из стека.
		lb = lbstack[stackpos];
		ub = ubstack[stackpos];
		stackpos--;

		do {
			// Шаг 1. Разделение по элементу pivot
			ppos = (lb + ub) >> 1;
			i = lb; j = ub; pivot = a[ppos];
			do {
				while (a[i] < pivot) i++;
				while (pivot < a[j]) j--;
				if (i <= j) {
					temp = a[i]; a[i] = a[j]; a[j] = temp;
					i++; j--;
				}
			} while (i <= j);

			// Сейчас указатель i указывает на начало правого подмассива,
			// j - на конец левого (см. иллюстрацию выше), lb ? j ? i ? ub.
			// Возможен случай, когда указатель i или j выходит за границу массива

			// Шаги 2, 3. Отправляем большую часть в стек и двигаем lb,ub
			if (i < ppos) { // правая часть больше
				if (i < ub) { // если в ней больше 1 элемента - нужно
					stackpos++; // сортировать, запрос в стек
					lbstack[stackpos] = i;
					ubstack[stackpos] = ub;
				}
				ub = j; // следующая итерация разделения
						// будет работать с левой частью
			}
			else { // левая часть больше
				if (j > lb) {
					stackpos++;
					lbstack[stackpos] = lb;
					ubstack[stackpos] = j;
				}
				lb = i;
			}
		} while (lb < ub); // пока в меньшей части более 1 элемента
	} while (stackpos != 0); // пока есть запросы в стеке
}

/// sort centers shapes by x axis
void qsort_centers_x(struct center_st2 *a, indexer size) {

	long long i, j; // указатели, участвующие в разделении
	long long lb, ub; // границы сортируемого в цикле фрагмента

	long long lbstack[MAX_STACK], ubstack[MAX_STACK]; // стек запросов
													// каждый запрос задается парой значений,
													// а именно: левой(lbstack) и правой(ubstack)
													// границами промежутка
	unsigned char stackpos = 1; // текущая позиция стека
	long long ppos; // середина массива
	struct center_st2 pivot; // опорный элемент
	struct center_st2 temp;

	lbstack[1] = 0;
	ubstack[1] = size - 1;

	//FILE *f1;
	//fopen_s(&f1, "c:/projects/tmp/1/cmp1.txt", "w");
	//unsigned __int64 clock1;
	//clock1 = __rdtsc();

	do {
		// Взять границы lb и ub текущего массива из стека.
		lb = lbstack[stackpos];
		ub = ubstack[stackpos];
		stackpos--;

		do {
			// Шаг 1. Разделение по элементу pivot
			ppos = (lb + ub) >> 1;
			i = lb;
			j = ub;
			memcpy(&pivot, &(a[ppos]), sizeof(struct center_st2));
			do {
				while (a[i].cx < pivot.cx) {
					//fprintf(f1, "i = %u, p = %u, diff = %f (x = %f, y = %f)\n", i, ppos, pivot.cx - a[i].cx, a[i].cx, a[i].cy);
					i++;
				}
				while (pivot.cx < a[j].cx) {
					//fprintf(f1, "JJ = %u, p = %u, diff = %f (x = %f, y = %f)\n", j, ppos, a[j].cx - pivot.cx, a[i].cx, a[i].cy);
					j--;
				}
				if (i <= j) {
					//temp = a[i]; a[i] = a[j]; a[j] = temp;
					if (i != j) {
						memcpy(&temp, &(a[i]), sizeof(struct center_st2));
						memcpy(&(a[i]), &(a[j]), sizeof(struct center_st2));
						memcpy(&(a[j]), &temp, sizeof(struct center_st2));
						//fprintf(f1, "change i = %u, j = %u, center1 = %f, center2 = %f\n", i, j, a[i].cx, a[j].cx);
					}
					i++; j--;
				}
			} while (i <= j);

			// Сейчас указатель i указывает на начало правого подмассива,
			// j - на конец левого (см. иллюстрацию выше), lb ? j ? i ? ub.
			// Возможен случай, когда указатель i или j выходит за границу массива

			// Шаги 2, 3. Отправляем большую часть в стек и двигаем lb,ub
			if (i < ppos) { // правая часть больше
				if (i < ub) { // если в ней больше 1 элемента - нужно
					stackpos++; // сортировать, запрос в стек
					lbstack[stackpos] = i;
					ubstack[stackpos] = ub;
				}
				ub = j; // следующая итерация разделения
						// будет работать с левой частью
			}
			else { // левая часть больше
				if (j > lb) {
					stackpos++;
					lbstack[stackpos] = lb;
					ubstack[stackpos] = j;
				}
				lb = i;
			}
		} while (lb < ub); // пока в меньшей части более 1 элемента
	} while (stackpos != 0); // пока есть запросы в стеке
							 //fclose(f1);

							 /*	char ch[64];
							 sprintf_s(ch, 64, "time = %lld", __rdtsc() - clock1);
							 lprintf(ch);
							 */
}

/// sort centers shapes by y axis
void qsort_centers_y(struct center_st2 *a, indexer size) {

	long long i, j; // указатели, участвующие в разделении
	long long lb, ub; // границы сортируемого в цикле фрагмента

	long long lbstack[MAX_STACK], ubstack[MAX_STACK]; // стек запросов
													// каждый запрос задается парой значений,
													// а именно: левой(lbstack) и правой(ubstack)
													// границами промежутка
	unsigned char stackpos = 1; // текущая позиция стека
	long long ppos; // середина массива
	struct center_st2 pivot; // опорный элемент
	struct center_st2 temp;

	lbstack[1] = 0;
	ubstack[1] = size - 1;

	//FILE *f1;
	//fopen_s(&f1, "c:/projects/tmp/1/cmp1.txt", "w");
	//unsigned __int64 clock1;
	//clock1 = __rdtsc();

	do {
		// Взять границы lb и ub текущего массива из стека.
		lb = lbstack[stackpos];
		ub = ubstack[stackpos];
		stackpos--;

		do {
			// Шаг 1. Разделение по элементу pivot
			ppos = (lb + ub) >> 1;
			i = lb;
			j = ub;
			memcpy(&pivot, &(a[ppos]), sizeof(struct center_st2));
			do {
				while (a[i].cy < pivot.cy) {
					//fprintf(f1, "i = %u, p = %u, diff = %f (x = %f, y = %f)\n", i, ppos, pivot.cx - a[i].cx, a[i].cx, a[i].cy);
					i++;
				}
				while (pivot.cy < a[j].cy) {
					//fprintf(f1, "JJ = %u, p = %u, diff = %f (x = %f, y = %f)\n", j, ppos, a[j].cx - pivot.cx, a[i].cx, a[i].cy);
					j--;
				}
				if (i <= j) {
					//temp = a[i]; a[i] = a[j]; a[j] = temp;
					if (i != j) {
						memcpy(&temp, &(a[i]), sizeof(struct center_st2));
						memcpy(&(a[i]), &(a[j]), sizeof(struct center_st2));
						memcpy(&(a[j]), &temp, sizeof(struct center_st2));
						//fprintf(f1, "change i = %u, j = %u, center1 = %f, center2 = %f\n", i, j, a[i].cx, a[j].cx);
					}
					i++; j--;
				}
			} while (i <= j);

			// Сейчас указатель i указывает на начало правого подмассива,
			// j - на конец левого (см. иллюстрацию выше), lb ? j ? i ? ub.
			// Возможен случай, когда указатель i или j выходит за границу массива

			// Шаги 2, 3. Отправляем большую часть в стек и двигаем lb,ub
			if (i < ppos) { // правая часть больше
				if (i < ub) { // если в ней больше 1 элемента - нужно
					stackpos++; // сортировать, запрос в стек
					lbstack[stackpos] = i;
					ubstack[stackpos] = ub;
				}
				ub = j; // следующая итерация разделения
						// будет работать с левой частью
			}
			else { // левая часть больше
				if (j > lb) {
					stackpos++;
					lbstack[stackpos] = lb;
					ubstack[stackpos] = j;
				}
				lb = i;
			}
		} while (lb < ub); // пока в меньшей части более 1 элемента
	} while (stackpos != 0); // пока есть запросы в стеке
							 //fclose(f1);

							 /*	char ch[64];
							 sprintf_s(ch, 64, "time = %lld", __rdtsc() - clock1);
							 lprintf(ch);
							 */
}

/// sort centers shapes by x axis
void qsort_node_centers_x(struct center_node_st *a, indexer size) {

	long long i, j; // указатели, участвующие в разделении
	long long lb, ub; // границы сортируемого в цикле фрагмента

	long long lbstack[MAX_STACK], ubstack[MAX_STACK]; // стек запросов
													// каждый запрос задается парой значений,
													// а именно: левой(lbstack) и правой(ubstack)
													// границами промежутка
	unsigned char stackpos = 1; // текущая позиция стека
	long long ppos; // середина массива
	struct center_node_st pivot; // опорный элемент
	struct center_node_st temp;

	lbstack[1] = 0;
	ubstack[1] = size - 1;

	//FILE *f1;
	//fopen_s(&f1, "c:/projects/tmp/1/cmp1.txt", "w");
	//unsigned __int64 clock1;
	//clock1 = __rdtsc();

	do {
		// Взять границы lb и ub текущего массива из стека.
		lb = lbstack[stackpos];
		ub = ubstack[stackpos];
		stackpos--;

		do {
			// Шаг 1. Разделение по элементу pivot
			ppos = (lb + ub) >> 1;
			i = lb;
			j = ub;
			memcpy(&pivot, &(a[ppos]), sizeof(struct center_node_st));
			do {
				while (a[i].cx < pivot.cx) {
					//fprintf(f1, "i = %u, p = %u, diff = %f (x = %f, y = %f)\n", i, ppos, pivot.cx - a[i].cx, a[i].cx, a[i].cy);
					i++;
				}
				while (pivot.cx < a[j].cx) {
					//fprintf(f1, "JJ = %u, p = %u, diff = %f (x = %f, y = %f)\n", j, ppos, a[j].cx - pivot.cx, a[i].cx, a[i].cy);
					j--;
				}
				if (i <= j) {
					//temp = a[i]; a[i] = a[j]; a[j] = temp;
					if (i != j) {
						memcpy(&temp, &(a[i]), sizeof(struct center_node_st));
						memcpy(&(a[i]), &(a[j]), sizeof(struct center_node_st));
						memcpy(&(a[j]), &temp, sizeof(struct center_node_st));
						//fprintf(f1, "change i = %u, j = %u, center1 = %f, center2 = %f\n", i, j, a[i].cx, a[j].cx);
					}
					i++; j--;
				}
			} while (i <= j);

			// Сейчас указатель i указывает на начало правого подмассива,
			// j - на конец левого (см. иллюстрацию выше), lb ? j ? i ? ub.
			// Возможен случай, когда указатель i или j выходит за границу массива

			// Шаги 2, 3. Отправляем большую часть в стек и двигаем lb,ub
			if (i < ppos) { // правая часть больше
				if (i < ub) { // если в ней больше 1 элемента - нужно
					stackpos++; // сортировать, запрос в стек
					lbstack[stackpos] = i;
					ubstack[stackpos] = ub;
				}
				ub = j; // следующая итерация разделения
						// будет работать с левой частью
			}
			else { // левая часть больше
				if (j > lb) {
					stackpos++;
					lbstack[stackpos] = lb;
					ubstack[stackpos] = j;
				}
				lb = i;
			}
		} while (lb < ub); // пока в меньшей части более 1 элемента
	} while (stackpos != 0); // пока есть запросы в стеке
							 //fclose(f1);

							 /*	char ch[64];
							 sprintf_s(ch, 64, "time = %lld", __rdtsc() - clock1);
							 lprintf(ch);
							 */
}

/// sort centers shapes by y axis
void qsort_node_centers_y(struct center_node_st *a, indexer size) {

	long long i, j; // указатели, участвующие в разделении
	long long lb, ub; // границы сортируемого в цикле фрагмента

	long long lbstack[MAX_STACK], ubstack[MAX_STACK]; // стек запросов
													// каждый запрос задается парой значений,
													// а именно: левой(lbstack) и правой(ubstack)
													// границами промежутка
	unsigned char stackpos = 1; // текущая позиция стека
	long long ppos; // середина массива
	struct center_node_st pivot; // опорный элемент
	struct center_node_st temp;

	lbstack[1] = 0;
	ubstack[1] = size - 1;

	//FILE *f1;
	//fopen_s(&f1, "c:/projects/tmp/1/cmp1.txt", "w");
	//unsigned __int64 clock1;
	//clock1 = __rdtsc();

	do {
		// Взять границы lb и ub текущего массива из стека.
		lb = lbstack[stackpos];
		ub = ubstack[stackpos];
		stackpos--;

		do {
			// Шаг 1. Разделение по элементу pivot
			ppos = (lb + ub) >> 1;
			i = lb;
			j = ub;
			memcpy(&pivot, &(a[ppos]), sizeof(struct center_node_st));
			do {
				while (a[i].cy < pivot.cy) {
					//fprintf(f1, "i = %u, p = %u, diff = %f (x = %f, y = %f)\n", i, ppos, pivot.cx - a[i].cx, a[i].cx, a[i].cy);
					i++;
				}
				while (pivot.cy < a[j].cy) {
					//fprintf(f1, "JJ = %u, p = %u, diff = %f (x = %f, y = %f)\n", j, ppos, a[j].cx - pivot.cx, a[i].cx, a[i].cy);
					j--;
				}
				if (i <= j) {
					//temp = a[i]; a[i] = a[j]; a[j] = temp;
					if (i != j) {
						memcpy(&temp, &(a[i]), sizeof(struct center_node_st));
						memcpy(&(a[i]), &(a[j]), sizeof(struct center_node_st));
						memcpy(&(a[j]), &temp, sizeof(struct center_node_st));
						//fprintf(f1, "change i = %u, j = %u, center1 = %f, center2 = %f\n", i, j, a[i].cx, a[j].cx);
					}
					i++; j--;
				}
			} while (i <= j);

			// Сейчас указатель i указывает на начало правого подмассива,
			// j - на конец левого (см. иллюстрацию выше), lb ? j ? i ? ub.
			// Возможен случай, когда указатель i или j выходит за границу массива

			// Шаги 2, 3. Отправляем большую часть в стек и двигаем lb,ub
			if (i < ppos) { // правая часть больше
				if (i < ub) { // если в ней больше 1 элемента - нужно
					stackpos++; // сортировать, запрос в стек
					lbstack[stackpos] = i;
					ubstack[stackpos] = ub;
				}
				ub = j; // следующая итерация разделения
						// будет работать с левой частью
			}
			else { // левая часть больше
				if (j > lb) {
					stackpos++;
					lbstack[stackpos] = lb;
					ubstack[stackpos] = j;
				}
				lb = i;
			}
		} while (lb < ub); // пока в меньшей части более 1 элемента
	} while (stackpos != 0); // пока есть запросы в стеке
							 //fclose(f1);

							 /*	char ch[64];
							 sprintf_s(ch, 64, "time = %lld", __rdtsc() - clock1);
							 lprintf(ch);
							 */
}


#endif // SORT_H_HEADERS