

/* Program "lookat": utility to examine subsequences of a FASTA file. If the
   start position is greater than the end position, the reverse complement of
   the sequence is given (obviously, this feature is only meaningful for 
   DNA/RNA sequences!). The "+" and "-" options allow the user to easily add
   downstream and upstream sequence respectively to the displayed sequence. 
   This is quite a useful feature for those of us bad at arithmetic! 

   Sample runs:
   lookat EC.contigs 1 40

   >ECOLI_1_40
   AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTG

   lookat EC.contigs 40 1

   >ECOLI_40_1
   CACAGAGACATATTGCCCGTTGCAGTCAGAATGAAAAGCT

   lookat +5 EC.contigs 1 40

   >ECOLI_1_45
   AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGAT 
   
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
  Sequence sequence, subseq;
  FILE *infile = NULL;
  int i, upstream = 0, downstream = 0, start = 0, end = 0;
  char *arg;
  char contigname[500] = "\0'";

  if (argc < 4)			/* need at least FASTA file, start, and end positions */
    usage ();

  /* process arguments */

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (arg[0] == '-')
      upstream = atoi (arg);
    else if (arg[0] == '+')
      downstream = atoi (arg);
    else if (infile == NULL) {
      infile = fopen (arg, "r");
      if (!infile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
	exit (1);
      }
    }
    else if ((contigname[0] == '\0') && (atoi (arg) == 0))
      strcpy (contigname, arg);
    else if (start == 0)
      start = atoi (arg);
    else if (end == 0)
      end = atoi (arg);
  }

  /* if contig-name specified, load that contig */
  if (contigname[0] != '\0')
    sequence = loadSpecificContig (contigname, infile);
  else				/* just load the first contig by default */
    sequence = loadNextContig (infile);
  if (sequence == NULL) {	/* oops -- contig not found */
    fprintf (stderr, "Cannot find contig %s\n", contigname);
    exit (1);
  }
  if (start < end)		/* if forward strand, upstream & downstream adj to be added */
    subseq = lookAt (sequence, start + upstream, end + downstream);
  else				/* subtract upstream & downstream adj */
    subseq = lookAt (sequence, start - upstream, end - downstream);
  /* now print out the subsequence */
  printFasta (subseq, stdout);
  /* free memory */
  freeSequence (sequence);
  freeSequence (subseq);
  exit (0);
}

/* display usage information */

void
usage (void)
{
  fprintf (stderr, "Usage: lookat [options] fasta-file [contig-name] start end\n");
  fprintf (stderr, "valid lookat options:\n");
  fprintf (stderr, "     -num include num bases upstream of given region\n");
  fprintf (stderr, "     +num include num bases downstream of given region\n");
  exit (1);
}
