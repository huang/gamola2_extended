/* 
    scanblastpairs.c -- Program for creating a triplets file from a 
                        blastpairs file
 
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
#include "CodonList.h"
#include "Sequence.h"
#include "BlastIndex.h"


void usage (void);
void scanBlastPairs (Sequence s, BlastIndex index, FILE * blastpairs,
		     FILE * triplets, int *changes);

int
main (int argc, char *argv[])
{
  FILE *seqfile, *blastpairs, *triplets;
  Sequence s;
  BlastIndex index;
  int changes[4] =
  {0, 0, 0, 0};

  if (argc < 4)
    usage ();

  seqfile = fopen (argv[1], "r");
  blastpairs = fopen (argv[2], "r");
  triplets = fopen (argv[3], "w");
  if (seqfile == NULL) {
    fprintf (stderr, "%s: can't open %s\n", argv[0], argv[1]);
    exit (1);
  }
  if (blastpairs == NULL) {
    fprintf (stderr, "%s: can't open %s\n", argv[0], argv[2]);
    exit (1);
  }
  if (triplets == NULL) {
    fprintf (stderr, "%s: can't open %s\n", argv[0], argv[3]);
    exit (1);
  }
  index = newBlastIndex (blastpairs);
  while (!feof (seqfile)) {
    s = loadNextContig (seqfile);
    scanBlastPairs (s, index, blastpairs, triplets, changes);
    freeSequence (s);
  }
  fprintf (triplets, "changes: %d %d %d %d\n", changes[0], changes[1],
	   changes[2], changes[3]);
  freeBlastIndex (index);
  fclose (seqfile);
  fclose (blastpairs);
  fclose (triplets);
  exit (0);
}

void 
scanBlastPairs (Sequence s, BlastIndex index, FILE * blastpairs,
		FILE * triplets, int *changes)
{
  BlastIndex currentBlast = index;
  Sequence seq1, seq2, antiseq;
  CodonList list = newCodonList (s->header, s->length);
  int i, diff, start, end, temp;

  while (currentBlast != NULL) {
    if (strcmp (s->header, currentBlast->name) == 0) {
      fseek (blastpairs, currentBlast->matchPos, SEEK_SET);
      start = currentBlast->start;
      end = currentBlast->end;
      seq1 = loadNextContig (blastpairs);
      seq2 = loadNextContig (blastpairs);
      if (start > end) {
	antiseq = complement (seq1);
	freeSequence (seq1);
	seq1 = antiseq;
	antiseq = complement (seq2);
	freeSequence (seq2);
	seq2 = antiseq;
	temp = start;
	start = end;
	end = temp;
      }
      if (seq1->length != seq2->length) {
	fprintf (stderr, "skipping %s %d-%d inconsistent length\n",
		 currentBlast->name, currentBlast->start, currentBlast->end);
      }
      else {
	for (i = 0; i < seq2->length - 2; i++) {
	  diff = codonDiff (codonAt (seq1, i + 1), codonAt (seq2, i + 1));
	  changes[diff]++;
	  if (diff > 0)
	    addCodon (list, codonAt (seq2, i + 1), start + i);
	}
      }
      freeSequence (seq1);
      freeSequence (seq2);
    }
    currentBlast = nextBlastNode (currentBlast);
  }
  printCodonList (list, triplets);
  freeCodonList (list);
}

void
usage (void)
{
  fprintf (stderr, "Usage: scanblastpairs seq blast-pairs triplets\n");
  exit (1);
}
