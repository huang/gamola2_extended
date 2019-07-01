/* 
    CodonList.c -- Data structure for compactly storing occurrences of
    various codons at a given sequence position 
 
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
#include <string.h>
#include "Sequence.h"
#include "CodonList.h"

CodonList 
newCodonList (char *name, int length)
{
  CodonList list= (CodonList) CodonListCalloc (1, sizeof (CodonListStruct));
  int i;

  list->list = (CodonNode **) CodonListCalloc (length, sizeof (CodonList *));
  list->length = length;
  list->name = strdup (name);
  for (i = 0; i < list->length; i++) {
    list->list[i] = NULL;
  }
  return list;
}

void 
addCodon (CodonList list, char *codon, int pos)
{
  CodonNode *currentNode = NULL, *prevNode = NULL;
  int cn = codonNumber (codon);

  if (cn < 0)
    return;
  pos--;
  if ((pos >= list->length) || (pos < 0))
    return;

  if (list->list[pos] == NULL) {
    list->list[pos] = (CodonNode *) CodonListCalloc (1, sizeof (CodonNode));
    (list->list[pos])->codon = (unsigned char) cn;
    (list->list[pos])->num = 1;
    (list->list[pos])->next = NULL;
    return;
  }

  currentNode = list->list[pos];
  while (currentNode != NULL) {
    if (cn == currentNode->codon) {
      currentNode->num++;
      return;
    }
    prevNode = currentNode;
    currentNode = currentNode->next;
  }
  prevNode->next = (CodonNode *) CodonListCalloc (1, sizeof (CodonNode));
  (prevNode->next)->codon = (unsigned char) cn;
  (prevNode->next)->num = 1;
  (prevNode->next)->next = NULL;
}

void 
printCodonList (CodonList list, FILE * file)
{
  int i;
  char *codon;
  CodonNode *currentNode;

  fprintf (file, "%s\n", list->name);
  for (i = 0; i < list->length; i++) {
    currentNode = list->list[i];
    if (currentNode != NULL) {
      fprintf (file, "%-7d", i + 1);
      while (currentNode != NULL) {
	codon = codonString (currentNode->codon);
	fprintf (file, " %3s %3d", codon, currentNode->num);
	free (codon);
	currentNode = currentNode->next;
      }
      fprintf (file, "\n");
    }
  }
}

void 
freeCodonList (CodonList list)
{
  int i;
  CodonNode *currentNode, *prevNode;

  free (list->name);
  for (i = 0; i < list->length; i++) {
    currentNode = list->list[i];
    while (currentNode != NULL) {
      prevNode = currentNode;
      currentNode = currentNode->next;
      free (prevNode);
    }
  }
  free(list);
}

void *
CodonListCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for CodonList\n");
    exit (1);
  }
  else
    return ptr;
}
