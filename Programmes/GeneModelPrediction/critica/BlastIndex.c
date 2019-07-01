/* 
    BlastIndex.c -- Data structure for indexing "blastpairs" files
 
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
#include <string.h>
#include <malloc.h>
#include "BlastIndex.h"


BlastIndex 
newBlastIndex (FILE * file)
{
  char line[500], locus[500];
  int start, end, counter = 0;
  BlastIndex newNode, curPtr, startPtr;
  double p;
  long pos;

  curPtr = startPtr = NULL;
  while (!feof (file)) {
    pos = ftell (file);
    fscanf (file, "%499[^\n]\n", line);
    if (line[0] != '>')
      continue;
    else
      counter++;
    if (counter % 2 == 1) {
      sscanf (line, ">%s %d-%d %*s %lf", locus, &start, &end, &p);
      newNode = (BlastIndex) BlastIndexCalloc ((size_t) 1, sizeof (BlastIndexNode));
      newNode->name = (char *)
	BlastIndexCalloc ((size_t) strlen (locus) + 1, sizeof (char));
      strcpy (newNode->name, locus);
      newNode->start = start;
      newNode->end = end;
      newNode->p = p;
      newNode->matchPos = pos;
      newNode->next = NULL;
      if (curPtr == NULL)
	startPtr = curPtr = newNode;
      else
	curPtr->next = newNode;
      curPtr = newNode;
    }
  }
  return startPtr;
}

BlastIndex 
nextBlastNode (BlastIndex index)
{
  if (index == NULL)
    return NULL;
  else
    return index->next;
}

void
freeBlastIndex (BlastIndex index)
{
  BlastIndex current = index, next;

  while (current != NULL) {
    next = current->next;
    free (current->name);
    free (current);
    current = next;
  }
}


void *
BlastIndexCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for BlastIndex\n");
    exit (1);
  }
  else
    return ptr;
}
