/* 
    empirical-matrix.c -- Program to creates a scoring matrix with 
                          conservative changes from a blast-pairs 
                          file created from a BLASTN run of known 
                          coding regions in the organism of interest.
 
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
#include <math.h>
#include "Sequence.h"
#include "GeneticCode.h"

typedef struct {
  int coding[4096];
  int noncoding[4096];
} Counts;

void fillCounts (Sequence seq1, Sequence seq2, Counts * counts);
void usage (void);

int
main (int argc, char *argv[])
{
  Sequence sequence1, sequence2;
  GeneticCode g = newGeneticCode (11);
  Counts counts;
  int i, codon1, codon2, totcod = 0, totnon = 0;
  char *codon1str, *codon2str, aa1, aa2;
  double score, fc, fnc;
  FILE *infile = NULL;

  if (argc != 2)
    usage ();

  infile = fopen (argv[1], "r");
  if (!infile) {
    fprintf (stderr, "%s: can't open %s\n", argv[0], argv[1]);
    exit (1);
  }
  for (i = 0; i < 4096; i++) {
    counts.coding[i] = 0;
    counts.noncoding[i] = 0;
  }
  while (!feof (infile)) {
    sequence1 = loadNextContig (infile);
    sequence2 = loadNextContig (infile);
    fillCounts (sequence1, sequence2, &counts);
    freeSequence (sequence1);
    freeSequence (sequence2);
  }
  fclose (infile);
  for (i = 0; i < 4096; i++) {
    if (counts.coding[i] == 0)
      counts.coding[i] = 1;
    if (counts.noncoding[i] == 0)
      counts.noncoding[i] = 1;
    totcod += counts.coding[i];
    totnon += counts.noncoding[i];
  }
  for (i = 0; i < 4096; i++) {
    codon1 = i / 64;
    codon2 = i % 64;
    fc = (double) counts.coding[i] / (double) totcod;
    fnc = (double) counts.noncoding[i] / (double) totnon;
    score = log (fc / fnc);
    codon1str = codonString (codon1);
    codon2str = codonString (codon2);
    aa1 = translateCodon (codon1str, g);
    aa2 = translateCodon (codon2str, g);
    fprintf (stdout, "%s%s %4d %c%c %8.3e\n", codon1str, codon2str,
	     (int) (score / 0.015), aa1, aa2, fnc);
    free (codon1str);
    free (codon2str);
  }
  freeGeneticCode (g);
  exit (0);
}

void 
fillCounts (Sequence seq1, Sequence seq2, Counts * counts)
{
  Sequence antiseq1 = complement (seq1);
  Sequence antiseq2 = complement (seq2);
  int i, codon1, codon2, anticodon1, anticodon2;

  for (i = 1; i < seq1->length - 1; i++) {
    codon1 = codonNumber (codonAt (seq1, i));
    codon2 = codonNumber (codonAt (seq2, i));
    anticodon1 = codonNumber (codonAt (antiseq1, i));
    anticodon2 = codonNumber (codonAt (antiseq2, i));
    if ((anticodon1 >= 0) && (anticodon2 >= 0)) {
      counts->noncoding[64 * anticodon1 + anticodon2]++;
      counts->noncoding[64 * anticodon2 + anticodon1]++;
    }
    if ((codon1 >= 0) && (codon2 >= 0)) {
      if (i % 3 == 1) {
	counts->coding[64 * codon1 + codon2]++;
	counts->coding[64 * codon2 + codon1]++;
      }
      else {
	counts->noncoding[64 * anticodon1 + anticodon2]++;
	counts->noncoding[64 * anticodon2 + anticodon1]++;
      }
    }
  }
  freeSequence (antiseq1);
  freeSequence (antiseq2);
}

void
usage (void)
{
  fprintf (stderr, "Usage: empirical-matrix coding-blast-pairs\n");
  exit (1);
}
