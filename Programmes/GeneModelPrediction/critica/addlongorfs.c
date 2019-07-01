/* 
    addlongorfs.c -- Program to add long open reading frames to a set of
                     coding regions
 
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Array.h"
#include "Sequence.h"
#include "GeneticCode.h"
#include "SeqMap.h"

typedef struct Orf {
  char *name;
  int start;
  int end;
  int dir;
  int len;
  int seqLength;
  char initSeq[4];
} Orf;

void addLongOrfs (Array orfs, Sequence seq, GeneticCode g,
		  int minLen, SeqMap map, int dir);
void initOrf (Orf * orf, int pos, char *name, int seqLength, int dir);
int compareOrfs (const void *a, const void *b);
void printOrf (const void *a, FILE * file);
int orfLen (Orf * orf);

void usage (void);


int
main (int argc, char *argv[])
{
  FILE *seqfile = NULL, *cdsfile = NULL;
  int i, j, start, end, validStart, minLen = 75;
  Sequence seq, antiSeq;
  GeneticCode g = newGeneticCode (11);
  Array orfs;
  Orf *orf;
  SeqMap map;
  char *arg, locus[255];

  if (argc < 2)
    usage ();

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (arg[0] == '-') {
      if (strstr (arg, "-orf-aa-length=") != NULL) {
	arg = strtok (arg, "=");
	arg = strtok (NULL, "=");
	minLen = atoi (arg);
      }
      else if (strstr (arg, "-genetic-code=") != NULL) {
	arg = strtok (arg, "=");
	setCode (g, atoi (strtok (NULL, "=")));
      }
      else {
	fprintf (stderr, "%s: option %s invalid\n", argv[0], arg);
	exit (1);
      }
    }
    else if (cdsfile == NULL) {
      cdsfile = fopen (arg, "r");
      if (cdsfile == NULL) {
	fprintf (stderr, "%s: %s not found.\n", argv[0], arg);
	exit (1);
      }
    }
    else if (seqfile == NULL) {
      seqfile = fopen (arg, "r");
      if (seqfile == NULL) {
	fprintf (stderr, "%s: %s not found.\n", argv[0], arg);
	exit (1);
      }
    }
  }
  if ((seqfile == NULL) || (cdsfile == NULL))
    usage ();

  while (!feof (seqfile)) {
    orfs=newArray(5000);
    seq = loadNextContig (seqfile);
    antiSeq = complement (seq);
    map = newSeqMap (seq->length);
    rewind (cdsfile);
    while (!feof (cdsfile)) {
      fscanf (cdsfile, "%s %d %d%*[^\n]", locus, &start, &end);
      if (strcmp (seq->header, locus) == 0)
	fillSeqMap (map, start, end, 1);
    }
    addLongOrfs (orfs, seq, g, minLen, map, 1);
    addLongOrfs (orfs, antiSeq, g, minLen, map, 2);
    sortArray(orfs, compareOrfs);
    for(j=0; j<orfs->elements; j++) {
      orf = (Orf *) orfs->array[j];
      validStart = 1;
      for (i = orf->start; i < orf->end; i += 3) {
	if (mapPos (map, i-1, orf->dir) > 0) {
	  validStart = 0;
	}
	if ((!validStart) && (isInitCodon (&seq->seq[i], g))) {
	  orf->start = i;
	  validStart = 1;
	  strncpy (orf->initSeq, &seq->seq[i], 3);
	  orf->initSeq[3] = '\0';
	}
      }
      if ((orf->len > minLen) && (validStart)) {
	printOrf (orf, stdout);
	fillSeqMap (map, orf->start, orf->end, 1);
      }
    }
    freeArray(orfs, free);
    freeSequence (seq);
    freeSequence (antiSeq);
    freeSeqMap (map);
  }
  fclose (seqfile);
  fclose (cdsfile);
  exit (0);
}


void
addLongOrfs (Array orfs, Sequence seq, GeneticCode g, int minLen, 
             SeqMap map, int dir)
{
  Orf *orf;
  int i, j, validStart;
  char aa;

  orf = (Orf *) SeqMapCalloc (1, sizeof (Orf));

  for (i = 1; i <= 3; i++) {
    initOrf (orf, i, seq->header, seq->length, dir);
    validStart = 1;

    for (j = i; j < seq->length-2; j += 3) {
      aa = translateCodon (&seq->seq[j-1], g);
      if ((!validStart) && (isInitCodon (&seq->seq[j-1], g))) {
	orf->start = j;
	validStart = 1;
	strncpy (orf->initSeq, &seq->seq[j-1], 3);
	orf->initSeq[3] = '\0';
      }
      if (mapPos (map, j, dir) > 0)
	validStart = 0;
      if ((aa == '*') || ((seq->length - j) < 3)) {
	orf->end = j + 2;
	orf->len = orfLen (orf);
	if ((orf->len > minLen) && (validStart == 1)) {
	  addToArray (orfs, orf);
	  orf = (Orf *) SeqMapCalloc (1, sizeof (Orf));
	}
	initOrf (orf, j, seq->header, seq->length, dir);
	validStart = 0;
      }
    }
  }
}

void 
initOrf (Orf * orf, int pos, char *name, int seqLength, int dir)
{
  orf->start = pos;
  orf->end = pos + 2;
  orf->name = name;
  orf->dir = dir;
  orf->seqLength = seqLength;
  orf->initSeq[0] = '-';
  orf->initSeq[1] = '\0';
}

int
compareOrfs (const void *a, const void *b)
{
  const Orf *orf1 = *(Orf **) a, *orf2 = *(Orf **) b;
  if (orf1->len < orf2->len)
    return 1;
  else if (orf1->len > orf2->len)
    return -1;
  else
    return 0;
}

int
orfLen (Orf * orf)
{
  return (1 + abs (orf->start - orf->end)) / 3;
}

void 
printOrf (const void *a, FILE * file)
{
  const Orf *orf = (Orf *) a;
  int start, end;
  double p = 1e-4 / orf->len;

  if (orf->dir == 1) {
    start = orf->start;
    end = orf->end;
  }
  else {
    start = orf->seqLength - orf->start +1;
    end = orf->seqLength - orf->end +1;
  }
  fprintf (file, "%10s %8d %8d %9.2e %4d %6d %6d %4d %3s %4d %3d %-10s\n",
	   orf->name, start, end, p, 0, 0, 0, 0, orf->initSeq, 0, 0, "-");
}

void
usage (void)
{
  fprintf (stderr, "Usage: addlongorfs [options] cds-file seq-file\n");
  fprintf (stderr, "Valid addlongorfs options:\n");
  fprintf (stderr, "\t-orf-aa-length=length\n");
  fprintf (stderr, "\t-genetic-code=number\n");
  exit (1);
}
