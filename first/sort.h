#ifndef SORT_H_HEADERS
#define SORT_H_HEADERS

#include "pch.h"

#define MAXSTACK 2048		// максимальный размер стека
///// http://algolist.manual.ru/sort/quick_sort.php

void q_sort_mas(coord *a, long size, unsigned *pos_idx) { //qSortI

  long i, j;   			// указатели, участвующие в разделении

  long lb, ub;  		// границы сортируемого в цикле фрагмента

  long lbstack[MAXSTACK], ubstack[MAXSTACK]; // стек запросов
                        // каждый запрос задается парой значений,
                        // а именно: левой(lbstack) и правой(ubstack) 
                        // границами промежутка

  long stackpos = 1;   	// текущая позиция стека
  long ppos;            // середина массива
  coord pivot;              // опорный элемент
  coord temp;
  unsigned temp_int; 

  lbstack[1] = 0;
  ubstack[1] = size-1;

  do {

    // Взять границы lb и ub текущего массива из стека.

    lb = lbstack[ stackpos ];
    ub = ubstack[ stackpos ];
    stackpos--;

    do {
      // Шаг 1. Разделение по элементу pivot

      ppos = ( lb + ub ) >> 1;
      i = lb; j = ub; pivot = a[ppos];

      do {
        while ( a[i] < pivot ) i++;
        while ( pivot < a[j] ) j--;

        if ( i <= j ) {
          temp = a[i]; a[i] = a[j]; a[j] = temp;
          temp_int = pos_idx[i]; pos_idx[i] = pos_idx[j]; pos_idx[j] = temp_int;
          i++; j--;
        }
      } while ( i <= j );

      // Сейчас указатель i указывает на начало правого подмассива,
      // j - на конец левого (см. иллюстрацию выше), lb ? j ? i ? ub.
      // Возможен случай, когда указатель i или j выходит за границу массива

      // Шаги 2, 3. Отправляем большую часть в стек  и двигаем lb,ub

      if ( i < ppos ) {     // правая часть больше

        if ( i < ub ) {     //  если в ней больше 1 элемента - нужно 
          stackpos++;       //  сортировать, запрос в стек
          lbstack[ stackpos ] = i;
          ubstack[ stackpos ] = ub;
        }
        ub = j;             //  следующая итерация разделения
                            //  будет работать с левой частью

      } else {       	    // левая часть больше

        if ( j > lb ) { 
          stackpos++;
          lbstack[ stackpos ] = lb;
          ubstack[ stackpos ] = j;
        }
        lb = i;
      }

    } while ( lb < ub );        // пока в меньшей части более 1 элемента

  } while ( stackpos != 0 );    // пока есть запросы в стеке
}

/// BAD SORT
void quick_sort_r_mass(coord *a, long N, unsigned *pos_idx) { // quickSortR - http://algolist.manual.ru/sort/quick_sort.php
// На входе - массив a[], a[N] - его последний элемент.

  long i = 0, j = N-1;    // поставить указатели на исходные места
  coord temp, p;
  unsigned temp_int;

  p = a[ N>>1 ];    // центральный элемент

  // процедура разделения
  do {
    while ( a[i] < p ) i++;
    while ( a[j] > p ) j--;

    if (i <= j) {
      temp = a[i]; a[i] = a[j]; a[j] = temp;
      temp_int = pos_idx[i]; pos_idx[i] = pos_idx[j]; pos_idx[j] = temp_int;
      i++; j--;
    }
  } while ( i<=j );


  // рекурсивные вызовы, если есть, что сортировать 
  if (j > 0) quick_sort_r_mass(a, j, pos_idx);
  if (N > i) quick_sort_r_mass(a + i, N - i, pos_idx);
}

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
// only for 8-bits digit
inline void min_max_sort(int *a, int *b) // https://habrahabr.ru/post/204682/
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

#endif // _WIN

#endif // SORT_H_HEADERS