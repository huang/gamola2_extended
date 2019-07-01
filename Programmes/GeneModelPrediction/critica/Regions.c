/* 
    Regions.c -- Routines for finding and manipulating regions of a sequence
                 enriched in coding support
 
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
#include "Array.h"
#include "Sequence.h"
#include "GeneticCode.h"
#include "ScoringMatrix.h"
#include "DicodonScores.h"
#include "Statistics.h"
#include "Scores.h"
#include "Regions.h"
#include "Features.h"

Array
findRegions (Sequence seq, Sequence antiSeq, Scores scores,
	     DicodonScores dicod, Statistics stats, GeneticCode g)
{
  Array array = newArray (10000);	/* arbitrary */

  findPeaks (array, seq, scores->compScores, scores->n, scores->nlabel,
	     dicod, stats, g, 1);
  findPeaks (array, antiSeq, scores->antiCompScores, scores->n, 
             scores->nlabel, dicod, stats, g, -1);

  return array;
}

Region *
initRegions (Sequence seq, Statistics stats, int dir, int n, int label)
{
  int i;
  Region *regions = (Region *) RegionCalloc (3, sizeof (Region));

  for (i = 0; i < 3; i++) {
    regions[i].name = seq->header;
    regions[i].start = -1;
    regions[i].seqLength = seq->length;
    regions[i].dir = (short) dir;
    regions[i].compScore = 0;
    regions[i].dicodonScore = 0;
    regions[i].initScore = 0;
    regions[i].sdScore = 0;
    regions[i].promScore = 0;
    regions[i].promErr = 9;
    regions[i].sdSeq[0] = '-';
    regions[i].sdSeq[1] = '\0';
    regions[i].initSeq[0] = '-';
    regions[i].initSeq[1] = '\0';
    regions[i].promSeq[0] = '-';
    regions[i].promSeq[1] = '\0';
    regions[i].sdPos = 0;
    regions[i].p = 1;
    regions[i].k = stats->k[n];
    regions[i].lambda = stats->lambda[n];
    regions[i].matrix = label;
    regions[i].n = n;
  }
  return regions;
}


void 
findPeaks (Array array, Sequence seq, short **compScores, int nval,
	   int *nlabel, DicodonScores dicod, Statistics stats,
	   GeneticCode g, int dir)
{
  int n, i;
  int score[3], max[3], combined[3];
  double dicodonScore[3];
  Region *regions;

  for (n = 0; n < nval; n++) {
    for (i = 0; i < 3; i++) {
      score[i] = 0;
      dicodonScore[i] = 0;
      max[i] = 0;
    }

    regions = initRegions (seq, stats, dir, n, nlabel[n]);
    for (i = 0; i < seq->length-2; i++) {
      score[i % 3] += compScores[n][i];
      if (i>=3)
        dicodonScore[i % 3] += scoreDicodon (&seq->seq[i-3], dicod);
      combined[i % 3] = score[i % 3] + (int) (dicodonScore[i % 3] / 
                                              stats->lambda[n]);

      /* score above zero, possible start of coding region */
      if ((combined[i % 3] > 0) && (regions[i % 3].start == -1))
	regions[i % 3].start = i;

      /* new high score */
      if (combined[i % 3] > max[i % 3]) {
	max[i % 3] = combined[i % 3];
	regions[i % 3].compScore = score[i % 3];
	regions[i % 3].dicodonScore = dicodonScore[i % 3];
	regions[i % 3].end = i + 2;
      }

      /* end of high-scoring segment */
      if ((combined[i % 3] < 0) || (isStopCodon (&seq->seq[i], g)) ||
	  (i >= seq->length - 3)) {
	regions[i % 3].p = computeP (regions[i % 3], stats);
	if (regions[i % 3].p < stats->threshold)
	  addToArray (array, copyRegion (regions[i % 3]));
	score[i % 3] = 0;
	dicodonScore[i % 3] = 0;
	max[i % 3] = 0;
	regions[i % 3].start = -1;
	regions[i % 3].compScore = 0;
	regions[i % 3].dicodonScore = 0;
      }
    }
    free (regions);
  }
}


/* returns the true (1-based, adjusted for forward or reverse strand) start of region */
int
realStart (Region * reg)
{
  if (reg->dir == 1) {
    return reg->start + 1;
  }
  else {
    return reg->seqLength - reg->start;
  }
}

/* returns the true (1-based, adjusted for forward or reverse strand) end of region */
int
realEnd (Region * reg)
{
  if (reg->dir == 1) {
    return reg->end + 1;
  }
  else {
    return reg->seqLength - reg->end;
  }
}


void
printRegion (Region * reg, Statistics stats, FILE * file)
{
  int dscore, iscore, sdscore, promscore;

  dscore = (int) (reg->dicodonScore / reg->lambda);
  iscore = (int) (reg->initScore / reg->lambda);
  sdscore = (int) (reg->sdScore / reg->lambda);
  promscore = (int) (reg->promScore / reg->lambda);

  if ((!stats->strict && (reg->p<0.99)) || (reg->p<stats->threshold)) {
    fprintf (file, "%10s %8d %8d %9.2e %4d %6d %6d %4d %3s", reg->name, 
             realStart (reg), realEnd (reg), reg->p, reg->matrix,
	     reg->compScore, dscore, iscore, reg->initSeq);
    fprintf(file," %4d %3d %10s", sdscore, reg->sdPos, reg->sdSeq);
    if (stats->promScan) {
      fprintf(file," %4d %3d %30s %2d", promscore, reg->promPos, 
              reg->promSeq, reg->promErr);
    }
    fprintf(file,"\n");
  }
}

void
extendRegions (Sequence seq, Sequence antiSeq, Array regions, Scores scores,
	       DicodonScores dicod, BonusScores initScores, 
               BonusScores sdScores, BonusScores promScores, 
               Statistics stats, GeneticCode g)
{
  int i;
  Region *reg;
  Sequence regSeq;
  short **compScores, **antiCompScores;
  char **matrix=dnaMatrix();

  sortArray (regions, sortRegionsByP);
  for (i = 0; i < regions->elements; i++) {
    reg = (Region *) regions->array[i];
    if (reg->dir == 1) {
      regSeq = seq;
      compScores = scores->compScores;
      antiCompScores = scores->antiCompScores;
    }
    else {
      regSeq = antiSeq;
      compScores = scores->antiCompScores;
      antiCompScores = scores->compScores;
    }
    *reg = findTerminator (*reg, regSeq, g);
    *reg = recomputeScore (*reg, compScores, regSeq, dicod, stats, g);
    findInitiator (*reg, compScores, antiCompScores, scores->n, regSeq,
		   dicod, initScores, sdScores, promScores, stats, g, matrix);
  }
  freeMatchingMatrix(matrix);
}

Region
recomputeScore (Region reg, short **compScores, Sequence seq,
		DicodonScores dicod, Statistics stats, GeneticCode g)
{
  Region newReg = reg;
  int i;

  newReg.compScore = 0;
  newReg.dicodonScore = 0;
  for (i = newReg.start; i < newReg.end; i += 3) {
    if (!isStopCodon (&seq->seq[i], g)) {
      newReg.dicodonScore += scoreDicodon (&seq->seq[i], dicod);
      newReg.compScore += compScores[reg.n][i];
      if (newReg.compScore + newReg.dicodonScore / newReg.lambda < 0) {
	newReg.compScore = 0;
	newReg.dicodonScore = 0;
      }
    }
  }
  newReg.p = computeP (newReg, stats);
  return newReg;
}

void
clearScores (Region reg, short **compScores, short **antiCompScores, int n)
{
  int i, j, sharedThirdPos;

  for (i = reg.start; i < reg.end - 2; i += 3) {
    for (j = 0; j < n; j++) {
      compScores[j][i] = 0;
      sharedThirdPos = reg.seqLength - i - 2;
      antiCompScores[j][sharedThirdPos] = 0;
    }
  }
}

void
findInitiator (Region reg, short **compScores, short **antiCompScores,
	       int n, Sequence seq, DicodonScores dicod, 
               BonusScores initScores, BonusScores sdScores, 
               BonusScores promScores, Statistics stats,
	       GeneticCode g, char **matrix)
{
  int i;
  Region newReg, bestReg;

  bestReg.p = 0.99;

  for (i = reg.start + 3; i < reg.end - 3; i += 3) {
    if (isInitCodon (&seq->seq[i], g)) {
      newReg = reg;
      newReg.start = i;
      initCopy(&newReg, seq, initScores);
      if (stats->sdScan)
        sdScan(&newReg, seq, matrix, sdScores);
      if (stats->promScan)
        promoterScan(&newReg, seq, matrix, promScores); 
      newReg = recomputeScore (newReg, compScores, seq, dicod, stats, g);
      if (newReg.p < bestReg.p)
	bestReg = newReg;
      if (newReg.p < bestReg.p * 100) 	
        printRegion (&newReg, stats, stdout);
    }
  }
  for (i = reg.start; i >= 0; i -= 3) {
    if (isStopCodon (&seq->seq[i], g))
      break;
    
    if ((isInitCodon (&seq->seq[i], g)) || (i<3)) {
      newReg = reg;
      newReg.start = i;
      if (isInitCodon (&seq->seq[i], g)) {
        initCopy(&newReg, seq, initScores);
        if (stats->sdScan)
          sdScan(&newReg, seq, matrix, sdScores);
        if (stats->promScan)
          promoterScan(&newReg, seq, matrix, promScores); 
      }
      else {
        offCopy(&newReg, initScores);
      }
      newReg = recomputeScore (newReg, compScores, seq, dicod, stats, g);
      if (newReg.p < bestReg.p) {
	bestReg.p = newReg.p;
	bestReg = newReg;
      }
      if (newReg.p < bestReg.p * 100)
	printRegion (&newReg, stats, stdout);
    }
  }
  if (bestReg.p < 0.99)
    clearScores (bestReg, compScores, antiCompScores, n);
}

Region
findTerminator (Region reg, Sequence seq, GeneticCode g)
{
  int i;

  for (i = reg.end - 2; i < reg.seqLength - 5; i += 3) {
    if (isStopCodon (&seq->seq[i], g))
      break;
    reg.end += 3;
  }
  return reg;
}

double
computeP (Region reg, Statistics stats)
{
  int score, sites = 2000;

  score = reg.compScore + (int) (((reg.dicodonScore * stats->alpha)
			    + reg.initScore + reg.sdScore + reg.promScore)
			   / reg.lambda);
  return -Nlm_Expm1 (-reg.k * sites * exp (-reg.lambda * score));
}



int
sortRegionsByP (const void *a, const void *b)
{
  const Region *reg1 = (Region *) a, *reg2 = (Region *) b;
  if (reg1->p > reg2->p)
    return 1;
  else if (reg1->p < reg2->p)
    return -1;
  else
    return 0;
}

int
sortRegionsByStart (const void *a, const void *b)
{
  Region *reg1 = (Region *) a, *reg2 = (Region *) b;
  int start1 = reg1->start;
  int start2 = reg2->start;

  if (start1 > start2)
    return 1;
  else if (start1 < start2)
    return -1;
  else
    return 0;
}


Region *
copyRegion (Region reg)
{
  Region *copy = (Region *) RegionCalloc (1, sizeof (Region));

  *copy = reg;
  return copy;
}

int
regionLen (Region reg)
{
  return (1 + abs (reg.start - reg.end)) / 3;
}

/* wrappers around Calloc to allow error checking */

void *
RegionCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for Region\n");
    exit (1);
  }
  else
    return ptr;
}
