/* Program "extractcoding" -- given a FASTA file of a genome and a list of 
   annotated genes, creates a FASTA file with a separate entry for each gene. 
   The output can include the gene descriptions if desired, and can also 
   translate the genes. 

   Sample runs:

   extractcoding EC.contigs EC.cds

   >ECOLI_190_255
   ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCG
   GGCTGA
   >ECOLI_337_2799
   ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTT
   GCCGATATTCTGGA ...

   extractcoding -prot EC.contigs EC.cds 

   >ECOLI_190_255 (translated)
   MKRISTTITTTITITTGNGAG*
   >ECOLI_337_2799 (translated)
   MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDA
   LPNISDAERIFAELL ...

   extractcoding -prot -desc EC.contigs EC.cds 

   >ECOLI_190_255 thrL b0001  thr operon leader peptide (translated)
   MKRISTTITTTITITTGNGAG*
   >ECOLI_337_2799 thrA b0002  aspartokinase I, homoserine dehydrogenase I (translated)
   MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDA
   LPNISDAERIFAELL ...

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
#include "GeneticCode.h"

void usage (void);

int
main (int argc, char *argv[])
{
  Sequence sequence, subseq, protseq;
  GeneticCode g = newGeneticCode (11);	/* bacterial genetic code by default */
  FILE *seqfile = NULL, *cdsfile = NULL;
  int i, start, end, outputProt = 0, addDesc = 0, min = 0, max = 0;
  char *arg, name[255], desc[255], newheader[255];

  if (argc < 3)			/* need at least FASTA file and cds-file */
    usage ();

  /* process arguments */

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (strstr (arg, "-prot") != NULL) {
      outputProt = 1;
    }
    else if (strstr (arg, "-desc") != NULL) {
      addDesc = 1;
    }
    else if (strstr (arg, "-min=") != NULL) {
      arg = strtok (arg, "=");
      min = atoi (strtok (NULL, "="));
    }
    else if (strstr (arg, "-max=") != NULL) {
      arg = strtok (arg, "=");
      max = atoi (strtok (NULL, "="));
    }
    else if (strstr (arg, "-genetic-code=") != NULL) {
      arg = strtok (arg, "=");
      setCode (g, atoi (strtok (NULL, "=")));
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

  /* load each sequence in turn */
  while (!feof (seqfile)) {
    sequence = loadNextContig (seqfile);
    rewind (cdsfile);
    while (!feof (cdsfile)) {	/* run through cds file */
      name[0] = 'X';
      name[1] = '\0';
      fscanf (cdsfile, "%s %d %d%[^\n]", name, &start, &end, desc);
      if (strcmp (sequence->header, name) == 0) {	/* are we in the right contig? */
	subseq = lookAt (sequence, start, end);		/* extract subsequence */
	if ((subseq->length / 3 < min) || ((max > 0)
					  && (subseq->length / 3 > max))) {
	  freeSequence (subseq);
	  continue;		/* too short or too long */
	}
	if (addDesc) {		/* should we add the description? */
	  strcpy (newheader, subseq->header);
	  strcat (newheader, desc);
	  setHeader (subseq, newheader);
	}
	if (outputProt) {	/* translate the gene? */
	  protseq = translate (subseq, g);
	  printFasta (protseq, stdout);
	  freeSequence (protseq);
	}
	else
	  printFasta (subseq, stdout);
	freeSequence (subseq);
      }
    }
    freeSequence (sequence);
  }
  fclose (seqfile);
  fclose (cdsfile);
  exit (0);
}

void
usage (void)
{
  fprintf (stderr, "Usage: extractcoding [options] fasta-file cds-file\n");
  fprintf (stderr, "valid extractcoding options:\n");
  fprintf (stderr, "     -desc  (include gene descriptions in headers)\n");
  fprintf (stderr, "     -prot  (output translated sequences)\n");
  fprintf (stderr, "     -genetic-code=NCBI-num\n");
  fprintf (stderr, "     -min=num (minimum length (in amino acids) of extracted genes)\n");
  fprintf (stderr, "     -max=num (maximum length (in amino acids) of extracted genes)\n");
  exit (1);
}
