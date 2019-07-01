/* 
    critica.c -- main program of CRITICA: Coding Region Identification Tool 
                Invoking Comparative Analysis
 
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
#include "Array.h"
#include "Sequence.h"
#include "GeneticCode.h"
#include "ScoringMatrix.h"
#include "DicodonScores.h"
#include "Statistics.h"
#include "Scores.h"
#include "Regions.h"
#include "Features.h"
#include "frameshifts.h"

void usage (void);

int
main (int argc, char *argv[])
{
  FILE *seqFile = NULL, *tripFile = NULL, *dicodonFile = NULL, *initFile,
   *sdFile = NULL, *promFile = NULL, *matrixFile = NULL;
  GeneticCode g = newGeneticCode (11);
  Array regions;
  Scores scores;
  DicodonScores dicodonScores;
  BonusScores initScores = NULL, sdScores = NULL, promScores = NULL;
  ScoringMatrix matrix;
  Statistics stats;
  Sequence seq, antiSeq;
  char *arg;
  int i, frameshiftThreshold = 0, quickStats=0, strict=0, nosd=0, prom=0;
  double alpha = -1, threshold = -1;

 
  if (argc < 3)
    usage ();

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (strstr (arg, "-genetic-code=") != NULL) {
      arg = strtok (arg, "=");
      setCode (g, atoi (strtok (NULL, "=")));
    }
    else if (strstr (arg, "-scoring-matrix=") != NULL) {
      arg = strtok (arg, "=");
      arg = strtok (NULL, "=");
      matrixFile = fopen (arg, "r");
      if (!matrixFile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
	exit (1);
      }
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
    else if (strstr (arg, "-init-scores=") != NULL) {
      arg = strtok (arg, "=");
      arg = strtok (NULL, "=");
      initFile = fopen (arg, "r");
      if (!initFile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
      }
      else {
	initScores = loadBonus (initFile);
	fclose (initFile);
      }
    }
    else if (strstr (arg, "-sd-scores=") != NULL) {
      arg = strtok (arg, "=");
      arg = strtok (NULL, "=");
      sdFile = fopen (arg, "r");
      if (!sdFile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
      }
      else {
	sdScores = loadBonus (sdFile);
	fclose (sdFile);
      }
    }
    else if (strstr (arg, "-prom-scores=") != NULL) {
      arg = strtok (arg, "=");
      arg = strtok (NULL, "=");
      promFile = fopen (arg, "r");
      if (!promFile) {
	fprintf (stderr, "%s: can't open %s\n", argv[0], arg);
      }
      else {
	promScores = loadBonus (promFile);
	fclose (promFile);
      }
    }
    else if (strstr (arg, "-frameshift-threshold=") != NULL) {
      arg = strtok (arg, "=");
      arg = strtok (NULL, "=");
      frameshiftThreshold = atoi (arg);
    }
    else if (strstr (arg, "-threshold=") != NULL) {
      arg = strtok (arg, "=");
      arg = strtok (NULL, "=");
      threshold = atof (arg);
    }
    else if (strstr (arg, "-strict-threshold") != NULL) {
      strict=1;
    }
    else if (strstr (arg, "-quick-stats") != NULL) {
      quickStats=1;
    }
    else if (strstr (arg, "-no-sdscores") != NULL) {
      nosd=1;
    }
    else if (strstr (arg, "-prom-find") != NULL) {
      prom=1;
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
    else {
      fprintf (stderr, "%s: argument %s not understood\n", argv[0], arg);
      exit (1);
    }
  }
  if (dicodonFile == NULL)	/* option -dicodon-scores not selected */
    dicodonScores = createEmptyDicodonScores ();
  else
    dicodonScores = loadDicodonScores (dicodonFile); /* load dicodon scores */
  if (matrixFile == NULL) 
    matrix = createDefaultScoringMatrix (g, tripFile);
  else
    matrix = loadMatrix(matrixFile);

  stats = newStatistics (matrix, dicodonScores, quickStats);

  if (alpha != -1)
    stats->alpha = alpha;
  if (threshold != -1)
    stats->threshold = threshold;
  if (strict == 1)
    stats->strict = 1;
  if (nosd == 1)
    stats->sdScan = 0;
  if (prom == 1)
    stats->promScan = 1;

  while (!feof (seqFile)) {
    seq = loadNextContig (seqFile);
    if (seq != NULL) { 
      antiSeq = complement(seq);
      scores = loadScores (tripFile, seq, antiSeq, matrix);
      regions = findRegions (seq, antiSeq, scores, dicodonScores, stats, g);
      if (frameshiftThreshold > 0) 
        frameshiftScan(regions, frameshiftThreshold);
      extendRegions (seq, antiSeq, regions, scores, dicodonScores, 
                     initScores, sdScores, promScores, stats, g);
      freeSequence (seq);
      freeSequence (antiSeq);
      freeScores (scores);
      freeArray (regions, free);
    }
  }
  fclose (seqFile);
  freeScoringMatrix (matrix);
  freeBonus (initScores);
  freeBonus (sdScores);
  freeBonus (promScores);
  freeDicodonScores (dicodonScores);
  freeStatistics (stats);
  exit (0);
}


void
usage (void)
{
  fprintf (stderr, "CRITICA 1.05 7/31/2000\n\n");
  fprintf (stderr, "Reference: Badger, Jonathan H. and Gary J. Olsen. CRITICA:\n");
  fprintf (stderr, "Coding Region Identification Tool Invoking Comparative Analysis.\n");
  fprintf (stderr, "Molecular Biology and Evolution 16: 512-524 (1999)\n\n");
  fprintf (stderr, "Usage: critica [options] seq triplet-scores\n");
  fprintf (stderr, "Valid critica options:\n");
  fprintf (stderr, "\t-scoring-matrix=file\n");
  fprintf (stderr, "\t-dicodon-scores=file\n");
  fprintf (stderr, "\t-init-scores=file\n");
  fprintf (stderr, "\t-sd-scores=file\n");
  fprintf (stderr, "\t-prom-scores=file\n");
  fprintf (stderr, "\t-no-sdscores\n");
  fprintf (stderr, "\t-prom-find\n");
  fprintf (stderr, "\t-genetic-code=number\n");
  fprintf (stderr, "\t-threshold=value\n");
  fprintf (stderr, "\t-alpha=value\n");
  fprintf (stderr, "\t-strict-threshold\n");
  fprintf (stderr, "\t-frameshift-threshold=value\n");
  fprintf (stderr, "\t-quick-stats\n");
  exit (1);
}
