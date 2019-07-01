/* 
    motiffind.c -- Program for finding a pattern in a FASTA file
 
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
#include "Sequence.h"

void usage (void);

int
main (int argc, char *argv[])
{
  Sequence sequence, subSeq, temp;
  FILE *infile = NULL;
  int i, pos, comp = 0, patLen = 0;
  char *arg, *pattern = NULL, **matrix=NULL;
  double gc;

  if (argc < 3)			/* need at least pattern and FASTA file */
    usage ();

  /* process arguments */

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (strstr (arg, "-compstrand") != NULL) {
      comp = 1;
    }
    else if (pattern == NULL) {
      pattern = strdup (arg);
      patLen = strlen (pattern);
    }
    else if (infile == NULL) {
      infile = fopen (arg, "r");
      if (!infile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
	exit (1);
      }
    }
  }
  if (infile == NULL)
    usage ();
  

  while (!feof (infile)) {
    sequence = loadNextContig (infile);
    if (comp) {
      temp = complement (sequence);
      freeSequence (sequence);
      sequence = temp;
    }
    if (matrix==NULL) {
      gc=gcContent(sequence);
      if (gc>0.2)
        matrix = dnaMatrix ();
      else
        matrix = identityMatrix ();
    }
    pos = findMotif (sequence, 1, pattern, matrix);
    while (pos > -1) {
      subSeq = lookAt (sequence, pos, pos + patLen - 1);
      if (!comp)
	fprintf (stdout, 
          "%s %d %s\n",sequence->header, pos, subSeq->seq);
      else
	fprintf (stdout, 
          "%s %d %s\n", sequence->header, 1 + sequence->length - pos,
	  subSeq->seq);
      freeSequence (subSeq);
      pos = findMotif (sequence, pos + 1, pattern, matrix);
    }
    /* free memory */
    freeSequence (sequence);
  }
  free (pattern);
  freeMatchingMatrix (matrix);
  exit (0);
}

/* display usage information */

void
usage (void)
{
  fprintf (stderr, "Usage: motiffind [-options] pattern fasta-file\n");
  fprintf (stderr, "Valid motiffind options:\n");
  fprintf (stderr, "\t-compstrand\n");
  exit (1);
}

