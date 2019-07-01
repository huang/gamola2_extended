/* Program "intergenic" -- finds all regions upstream of coding regions that
   don't overlap (in any frame, including the opposite strand), another
   coding region.

   sample runs:

   >ECOLI_1_189
   AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC
   TGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGG
   TCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC
   >ECOLI_256_336
   CGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACC
   AAAGGTAACGAGGTAACAACC ...

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
#include "SeqMap.h"

void usage (void);

int
main (int argc, char *argv[])
{
  Sequence sequence, subseq;
  FILE *seqfile = NULL, *cdsfile = NULL;
  int i, start, end, rstart, rend, min = 10, max = 0, addDesc = 0;
  char *arg, locus[255], desc[255], newheader[255];
  SeqMap map;

  if (argc < 3)			/* need at least FASTA file and cds-file */
    usage ();

  /* process arguments */

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (strstr (arg, "-min=") != NULL) {
      arg = strtok (arg, "=");
      min = atoi (strtok (NULL, "="));
    }
    else if (strstr (arg, "-max=") != NULL) {
      arg = strtok (arg, "=");
      max = atoi (strtok (NULL, "="));
    }
    else if (strstr (arg, "-desc") != NULL) {
      addDesc = 1;
    }
    else if (seqfile == NULL) {
      seqfile = fopen (arg, "r");
      if (!seqfile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
	exit (1);
      }
    }
    else if (cdsfile == NULL) {
      cdsfile = fopen (arg, "r");
      if (!cdsfile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
	exit (1);
      }
    }
  }

  if ((seqfile == NULL) || (cdsfile == NULL))
    usage ();

  while (!feof (seqfile)) {
    sequence = loadNextContig (seqfile);
    map = newSeqMap (sequence->length);
    rewind (cdsfile);  /* after map, created, need to step through each cds */
    while (!feof (cdsfile)) {
      fscanf (cdsfile, "%s %d %d\n", locus, &start, &end);
      if (strcmp (sequence->header, locus) == 0)	/* right contig */
	fillSeqMap (map, start, end, 1);
    }
    rewind (cdsfile);
    while (!feof (cdsfile)) {
      locus[0] = 'X';
      locus[1] = '\0';
      fscanf (cdsfile, "%s %d %d%[^\n]", locus, &start, &end, desc);
      if (strcmp (sequence->header, locus) == 0) {
	if (start < end) {
	  rend = start - 1;
	  rstart = rend;
	  while ((mapPos (map, rstart - 2, 1) == 0) && (rstart > 1))
	    rstart--;
	}
	else {
	  rend = start + 1;
	  rstart = rend;
	  while ((mapPos (map, rstart, 1) == 0) 
          && (rstart < sequence->length-1))
	    rstart++;
	}
	subseq = lookAt (sequence, rstart, rend);
	if ((subseq->length >= min) && ((subseq->length <= max) 
           || (max == 0))) {
	  if (addDesc) {	/* should we add the description? */
	    strcpy (newheader, subseq->header);
	    strcat (newheader, " ");
	    strcat (newheader, desc);
	    strcat (newheader, " (upstream region)");
	    setHeader (subseq, newheader);
	  }
	  printFasta (subseq, stdout);
	}
	freeSequence (subseq);
      }
    }
    freeSequence (sequence);
    freeSeqMap (map);
  }
  exit (0);
}


void
usage (void)
{
  fprintf (stderr, "Usage: intergenic [options] fasta-file cds-file\n");
  fprintf (stderr, "valid intergenic options:\n");
  fprintf (stderr, "     -desc  (include gene descriptions in headers)\n");
  fprintf (stderr, "     -min=num (minimum length (in amino acids) of extracted genes)\n");
  fprintf (stderr, "     -max=num (maximum length (in amino acids) of extracted genes)\n");
  exit (1);
}
