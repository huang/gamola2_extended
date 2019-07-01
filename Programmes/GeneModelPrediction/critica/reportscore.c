/* 
    reportscore.c -- Program for reporting the score of a potential coding
                     region 
 
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
#include <time.h>
#include "Sequence.h"
#include "GeneticCode.h"
#include "ScoringMatrix.h"
#include "DicodonScores.h"
#include "Scores.h"


void usage (void);

int
main (int argc, char *argv[])
{
  FILE *seqFile = NULL, *tripFile = NULL, *dicodonFile = NULL;
  GeneticCode g = newGeneticCode (11);
  Scores scores;
  Sequence seq, antiSeq;
  DicodonScores dicodonScores;
  ScoringMatrix matrix;
  char *arg, *contig=NULL;
  int i, start = 0 , end =0 ;
  double alpha = -1;
  double score=0;

  if (argc < 6)
    usage ();

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (strstr (arg, "-genetic-code=") != NULL) {
      arg = strtok (arg, "=");
      setCode (g, atoi (strtok (NULL, "=")));
    }
    else if (strstr (arg, "-dicodon-scores=") != NULL) {
      arg = strtok (arg, "=");
      arg = strtok (NULL, "=");
      dicodonFile = fopen (arg, "r");
      if (!dicodonFile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
	exit (1);
      }
    }
    else if (strstr (arg, "-alpha=") != NULL) {
      arg = strtok (arg, "=");
      arg = strtok (NULL, "=");
      alpha = atof (arg);
    }
    else if (seqFile == NULL) {
      seqFile = fopen (arg, "r");
      if (!seqFile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
	exit (1);
      }
    }
    else if (tripFile == NULL) {
      tripFile = fopen (arg, "r");
      if (!tripFile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
	exit (1);
      }
    }
    else if (contig==NULL) 
      contig=arg;
    else if (start==0) 
      start=atoi(arg);
    else if (end==0) 
      end=atoi(arg);
    else {
      fprintf (stderr, "%s: argument %s not understood\n", argv[0], arg);
      exit (1);
    }
  }
  if (dicodonFile == NULL)	/* option -dicodon-scores not selected */
    dicodonScores = createEmptyDicodonScores ();
  else
    dicodonScores = loadDicodonScores (dicodonFile); /* load dicodon scores */
  matrix = createDefaultScoringMatrix (g, tripFile);

  seq = loadSpecificContig (contig, seqFile);
  antiSeq = complement(seq);
  scores = loadScores (tripFile, seq, antiSeq, matrix);
  if (start< end) {
    for(i=start;i<end;i+=3) {
      score+=(scores->compScores[0][i-1]) +(scoreDicodon(&seq->seq[i-1], 
            dicodonScores)/0.015);

      if (score<0) score=0;
      printf("%7d %3d\n",i+1,(int)score);
    }
  }
  fclose (seqFile);
  freeSequence(seq);
  freeSequence(antiSeq);
  freeScores(scores);
  freeScoringMatrix (matrix);
  freeDicodonScores (dicodonScores);
  exit (0);
}


void
usage (void)
{
  fprintf (stderr, "Usage: reportscore [options] seq triplet-scores contig start end\n");
  fprintf (stderr, "Valid options:\n");
  fprintf (stderr, "\t-scoring-matrix=file\n");
  fprintf (stderr, "\t-dicodon-scores=file\n");
  fprintf (stderr, "\t-genetic-code=number\n");
  fprintf (stderr, "\t-alpha=value\n");
  exit (1);
}
