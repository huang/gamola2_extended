/* 
   Array.c -- dynamic array data structure for C

   Copyright (C) 1999  Jonathan H. Badger

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include "Array.h"

Array
newArray (int maxSize)
{
  Array array = (Array) malloc (sizeof (ArrayStruct));

  if (array == NULL) {
    fprintf (stderr, "not enough memory to allocate Array\n");
    exit (1);
  }
  else {
    array->maxSize = maxSize;
    array->elements = 0;
    array->dirty = 0;
    array->array = (void **) calloc (maxSize, sizeof (void *));
    if (array->array == NULL) {
      fprintf (stderr, "not enough memory to allocate Array\n");
      exit (1);
    }
  }
  return array;
}

void
addToArray (Array array, void *element)
{
  array->array[array->elements++] = element;
  array->dirty = 1;
  if (array->elements == array->maxSize) {
    array->maxSize *= 2;
    array->array = (void **) realloc (array->array, array->maxSize
				      * sizeof (void *));
    if (array->array == NULL) {
      fprintf (stderr, "not enough memory to allocate Array\n");
      exit (1);
    }
  }
}

void
freeArray (Array array, void (*freefunc) (void *))
{
  int i;

  for (i = 0; i < array->elements; i++) {
    freefunc (array->array[i]);
  }
  free (array->array);
  free (array);
}

/* merge routine for sortArray, below */
void
merge (void **array, void **work, int first, int middle, int last,
       int (*compare) (const void *, const void *))
{
  int i, j, k, n1;
  int n = last - first + 1;

  for (i = first, j = 0; i <= last;) {
    work[j++] = array[i++];
  }

  if (middle > last)
    middle = (first + last) / 2;

  n1 = middle - first + 1;

  for (i = first, j = 0, k = n1; i <= last; i++) {
    if ((j < n1) && ((k == n) || (compare (work[j], work[k]) < 0)))
      array[i] = work[j++];
    else
      array[i] = work[k++];
  }
}


void
sortArray (Array array, int (*compare) (const void *, const void *))
{
  /* iterative merge sort */

  void **work;
  int nt2 = array->elements * 2;
  int nm1 = array->elements - 1;
  int first, last, size;

  /* fprintf (stderr, "number of elements %d\n", array->elements); */

  work = (void **) calloc (array->elements, sizeof (void *));
  if (work == NULL && array->elements != 0 ) {
    fprintf (stderr, "not enough memory to allocate work array\n");
    exit (1);
  }
  for (size = 2; size < nt2; size *= 2) {
    for (first = 0; first < array->elements; first += size) {
      last = first + size - 1;
      merge (array->array, work, first, (first + last) / 2,
	     last < array->elements ? last : nm1, compare);
    }
  }
  free (work);
  array->dirty = 0;
}


/* sample code for Array routines */

/*
int *
copyInt (int i)
{
  int *copy = malloc (sizeof (int));

  *copy = i;
  return copy;
}

int
compare (const void *a, const void *b)
{
  const int int1 = **(int **) a, int2 = **(int **) b;
  if (int1 > int2)
    return 1;
  else if (int1 < int2)
    return -1;
  else
    return 0;
}

void
printElement (void *element, FILE * file)
{
  int i = *(int *) element;

  fprintf (file, "%d\n", i);
}

int
main ()
{
  Array array = newArray (10);
  int i;

  for (i = 20; i > 0; i--)
    addToArray (array, copyInt (i));

  for (i = 0; i < array->elements; i++) {
    printElement (array->array[i], stdout);
  }
  printf ("\n\n");
  sortArray (array, compare);
  for (i = 0; i < array->elements; i++) {
    printElement (array->array[i], stdout);
  }
  freeArray (array, free);
  exit (0);
}
*/
