/* 
    SeqMap.c -- Data structure for holding a map of a sequence to 
                allow marking of areas that are  coding. This is 
                needed for cases in which frequencies of codons 
                coding vs non-coding contexts need to be computed, 
                as well as finding intergenic regions 

 
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
#include "SeqMap.h"

/* creates a new map with the given length */

SeqMap 
newSeqMap (int length)
{
  SeqMap map= (SeqMap) SeqMapCalloc (length, sizeof (SeqMapStruct));
  map->coding = (char *) SeqMapCalloc (length, sizeof (char));
  map->length = length;
  return map;
}

/* frees memory used by SeqMap */

void 
freeSeqMap (SeqMap map)
{
  free (map->coding);
  free (map);
}

/* fill a SeqMap from start to end-> If AllFrames is 1, mark all frames
   in the region as coding-> Otherwise, just mark the translated frame */

void 
fillSeqMap (SeqMap map, int start, int end, int AllFrames)
{
  int i, step;

  if (AllFrames)
    step = 1;
  else
    step = 3;

  if (start > map->length-1)
    start = map->length-1;
  else if (start < 0)
    start = 0;
  if (end > map->length - 1)
    end = map->length - 1;
  else if (end < 0)
    end = 0;
  if (start < end)
    for (i = start-1; i < end; i += step)
      map->coding[i] = (char) (map->coding[i] | 1);
  else
    for (i = start-1; i >= end-1; i-= step)
      map->coding[i] = (char) (map->coding[i] | 2);
}

int 
mapPos (SeqMap map, int pos, int dir)
{
  if (dir == 1)
    return map->coding[pos];
  else
    return map->coding[map->length-pos-1];
}


/* wrapper around Calloc */

void *
SeqMapCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for SeqMap\n");
    exit (1);
  }
  else
    return ptr;
}
