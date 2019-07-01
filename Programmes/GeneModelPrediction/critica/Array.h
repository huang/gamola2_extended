#ifndef __ARRAY_H
#define __ARRAY_H

typedef struct {
  int elements;
  int maxSize;
  int dirty;
  void **array;
} ArrayStruct;

typedef ArrayStruct *Array;

Array newArray(int maxSize);
void addToArray(Array array, void * element);
void freeArray(Array array,void (*freefunc) (void *));
void sortArray(Array array, int (*compare)(const void *, const void *));

#endif /* __ARRAY_H */
