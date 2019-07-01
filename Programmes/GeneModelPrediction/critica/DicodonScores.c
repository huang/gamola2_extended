/* 
    DicodonScores.c -- Data structure for storing non-comparative
    dicodon scores
 
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
#include <malloc.h>
#include "Sequence.h"
#include "DicodonScores.h"

DicodonScores 
createEmptyDicodonScores ()
{
  DicodonScores dicod= (DicodonScores) DicodonScoresCalloc (1, 
                                       sizeof (DicodonScoresStruct));
  int i;
  double randProb = 1.0 / 4096;

  dicod->score = (double *) DicodonScoresCalloc (4097, sizeof (double));
  dicod->freq = (double *) DicodonScoresCalloc (4097, sizeof (double));
  for (i = 0; i < 4096; i++) {
    dicod->freq[i] = randProb;
  }
  dicod->score[4096]=0;
  dicod->minScore = 0;
  dicod->maxScore = 0;
  return dicod;
}

DicodonScores 
loadDicodonScores (FILE *file) 
{
  DicodonScores dicod= (DicodonScores) DicodonScoresCalloc (1, 
                                       sizeof (DicodonScoresStruct));
  char dicodonString[7];
  double score, freq;
  int cn1, cn2;


  dicod->minScore=1000;
  dicod->maxScore=-1000;
  dicod->score = (double *) DicodonScoresCalloc (4097, sizeof (double));
  dicod->freq = (double *) DicodonScoresCalloc (4097, sizeof (double));
  while(!feof(file)) {
    fscanf(file,"%s %le %le\n",dicodonString,&score,&freq);
    cn1=codonNumber(dicodonString);
    cn2=codonNumber(&dicodonString[3]);
    if ((cn1>=0) && (cn2>=0)) {
      dicod->score[(64*cn1)+cn2]=score;
      dicod->freq[(64*cn1)+cn2]=freq;
      if (score>dicod->maxScore) dicod->maxScore=score;
      if (score<dicod->minScore) dicod->minScore=score;
    }
  }
  return dicod;
}

void freeDicodonScores(DicodonScores dicod) {
  free(dicod->score);
  free(dicod->freq);
  free(dicod);
}

int dicodonNum(char *dicodon) {
  int cn1, cn2; 
  
  cn1=codonNumber(dicodon);
  cn2=codonNumber(&dicodon[3]);
  
  if ((cn1>=0) && (cn2>=0)) 
    return 64*cn1+cn2;
  else
    return 4096;
}

double scoreDicodon(char *dicodon, DicodonScores dicod) {
  return dicod->score[dicodonNum(dicodon)];
}

/* wrapper around Calloc to allow error checking */

void *
DicodonScoresCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for DicodonScores\n");
    exit (1);
  }
  else
    return ptr;
}
