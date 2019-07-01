/*
    dicodontable.c -- Program to calculate noncomparative dicodon scores

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
#include <math.h>
#include "Sequence.h"
#include "SeqMap.h"

typedef struct {
  double dicodon[64][64];
  double codon[64];
  double total_dicodons;
} coding_info;


void usage (void);
void count (SeqMap map, Sequence seq, coding_info * coding,
	    coding_info * noncoding, int dir);
void create_table (coding_info * coding, coding_info * noncoding,
		   double fraction);
void add_pseudocounts (coding_info * info);


int
main (int argc, char *argv[])
{
  FILE *seqfile = NULL, *cdsfile = NULL;
  int i, start, end;
  double fraction = 0;
  Sequence seq, antiSeq;
  SeqMap map;
  char *arg, locus[255];
  coding_info *coding, *noncoding;

  if (argc < 2)
    usage ();

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (arg[0] == '-') {
      if (strstr (arg, "-fraction-coding=") != NULL) {
	arg = strtok (arg, "=");
	arg = strtok (NULL, "=");
	fraction = atof (arg) / 6;
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

  coding    = SeqMapCalloc (1, sizeof (coding_info));
  noncoding = SeqMapCalloc (1, sizeof (coding_info));
  while (!feof (seqfile)) {
    seq = loadNextContig (seqfile);
    antiSeq = complement (seq);
    map = newSeqMap (seq->length);
    rewind (cdsfile);
    while (!feof (cdsfile)) {
      fscanf (cdsfile, "%s %d %d%*[^\n]\n", locus, &start, &end);
      if (strcmp (seq->header, locus) == 0) {
	fillSeqMap (map, start, end, 0);
      }
    }
    count (map, seq,     coding, noncoding, 1);
    count (map, antiSeq, coding, noncoding, 2);
    freeSequence (seq);
    freeSequence (antiSeq);
    freeSeqMap (map);
  }
  create_table (coding, noncoding, fraction);
  fclose (seqfile);
  fclose (cdsfile);
  free (coding);
  free (noncoding);
  exit (0);
}

void
count (SeqMap map, Sequence seq, coding_info * coding,
       coding_info * noncoding, int dir)
{
  int i, cn1, cn2;

  for (i = 0; i < seq->length - 5; i++) {
    cn1 = codonNumber (&seq->seq[i]);
    cn2 = codonNumber (&seq->seq[i+3]);
    if ((cn1 < 0) || (cn2 < 0)) continue;

    /* simplified by deferring the totals until the count is done -- GJO */
    /* mapPos returns 0 on noncoding
     *                1 on forward
     *                2 on reverse
     *                3 on both
     */

    if (mapPos(map, i, dir) & dir) coding->dicodon[cn1][cn2]++;
    else                           noncoding->dicodon[cn1][cn2]++;
  }
}

/*
 * Input:
 *
 *  fraction = prior probability of a hexamer coding
 *  coding->dicodon[i][j] = number of hexamers ij explicitly called coding
 *  noncoding->dicodon[i][j] = number of hexamers ij not classified (noncoding)
 *  coding->codon[i] = num. of expected coding hexamers that start with i
 *  noncoding->codon[i] = num. of expected noncoding hexamers that start with i
 *  coding->total_dicodons = total num. of dicodons explicitly called coding
 *  noncoding->total_dicodons = total num. of dicodons not classified
 *
 * Correct for unclassifed hexamers:
 *
 *  ratio = ratio of expected coding hexamers to called coding hexamers
 *  total_noncoding = total number of noncoding hexamers
 *  sum = coding->dicodon[i][j] + noncoding->dicodon[i][j]
 *  ratio = fraction * total_dicodons / total_called
 *
 *  coding->dicodon[i][j] = ratio * coding->dicodon[i][j]
 *  noncoding->dicodon[i][j] = sum - noncoding->dicodon[i][j]
 *
 */

void
create_table (coding_info * coding, coding_info * noncoding, double fraction)
{
  int i, j;
  char *codon1, *codon2;
  double ratio, deltaij, freq, score;
  double total_dicodons;

  /* do the initial totals -- GJO */

  coding->total_dicodons = 0;
  noncoding->total_dicodons = 0;
  for (i = 0; i < 64; i++) {
    for (j = 0; j < 64; j++) {
      coding->codon[i]    += coding->dicodon[i][j];
      noncoding->codon[i] += noncoding->dicodon[i][j];
    }
    coding->total_dicodons    += coding->codon[i];
    noncoding->total_dicodons += noncoding->codon[i];
  }

  total_dicodons = coding->total_dicodons + noncoding->total_dicodons;

  if (fraction > 0) {
    ratio = fraction * (total_dicodons / coding->total_dicodons);

    /* correct printing of initial fraction coding -- GJO */

    fprintf (stderr, "Initial fraction coding %5.3f. Extrapolating to %5.3f\n",
      6 * coding->total_dicodons / total_dicodons, 6 * fraction);
  }
  else
    ratio = 1.0;

  for (i = 0; i < 64; i++) {
    for (j = 0; j < 64; j++) {

      /* changed so that noncoding->dicodon cannot go negative -- GJO */
      /* defer new totals until all adjustments are all done -- GJO */

      deltaij = (ratio - 1) * coding->dicodon[i][j];
      if (deltaij > noncoding->dicodon[i][j]) deltaij = noncoding->dicodon[i][j];
      coding->dicodon[i][j]    += deltaij;
      noncoding->dicodon[i][j] -= deltaij;
    }
  }

  add_pseudocounts (coding);
  add_pseudocounts (noncoding);

  /* now rodo the totals -- GJO */

  coding->total_dicodons = 0;
  noncoding->total_dicodons = 0;
  for (i = 0; i < 64; i++) {
    coding->codon[i]    = 0;
    noncoding->codon[i] = 0;
    for (j = 0; j < 64; j++) {
      coding->codon[i]    += coding->dicodon[i][j];
      noncoding->codon[i] += noncoding->dicodon[i][j];
    }
    coding->total_dicodons    += coding->codon[i];
    noncoding->total_dicodons += noncoding->codon[i];
  }

  for (i = 0; i < 64; i++) {
    codon1 = codonString (i);
    for (j = 0; j < 64; j++) {
      codon2 = codonString (j);

      freq  = noncoding->dicodon[i][j] / noncoding->total_dicodons;
      score = log((coding->dicodon[i][j] / coding->codon[i]) /
               (noncoding->dicodon[i][j] / noncoding->codon[i]));

      /*  format changed from e to f for easier reading, more
       *  precision in the same space, and ease of sorting -- GJO
       */

      printf ("%3s%3s  %10.6f  %10.8f\n", codon1, codon2, score, freq);

      free (codon2);
    }
    free (codon1);
  }
}


#define PSEUDOCOUNT 0.5

/* if counts < PSEUDOCOUNT, counts=PSEUDOCOUNT
   (to avoid dividing by zero) */

void
add_pseudocounts (coding_info * info)
{
  int i, j;

  for (i = 0; i < 64; i++) {
    for (j = 0; j < 64; j++) {
      if (info->dicodon[i][j] < PSEUDOCOUNT) info->dicodon[i][j] = PSEUDOCOUNT;

      /* bookkeeping on totals is deferred -- GJO */
    }
  }
}

void
usage (void)
{
  fprintf (stderr, "Usage: dicodontable [options] cds-file seq-file\n");
  fprintf (stderr, "Valid dicontable options:\n");
  fprintf (stderr, "\t-fraction-coding=fraction\n");
  exit (1);
}
