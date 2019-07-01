/* 
    ScoringMatrix.c -- Data structure for holding a lookup table of
                       comparative scores for every pair of codons 
                       which CompScores uses to calculate the comparative 
                       score for each position
 
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
#include <malloc.h>
#include <ctype.h>
#include <string.h>
#include "Sequence.h"
#include "GeneticCode.h"
#include "ScoringMatrix.h"

/* creates the default set of scoring matrices as described in the CRITICA
   paper */


ScoringMatrix
createDefaultScoringMatrix (GeneticCode g, FILE * tripfile)
{
  ScoringMatrix matrix = (ScoringMatrix) ScoringMatrixCalloc (1, 
                                         sizeof (ScoringMatrixStruct));
  int n[5] =
  {8, 16, 32, 64, 128};
  int scores[] =
  {0, 0, 76, -99, 218, -48, 243, -28,	/* n=8 */
   0, 0, 65, -59, 197, -31, 230, -22,	/* n=16 */
   0, 0, 52, -36, 164, -16, 207, -14,	/* n=32 */
   0, 0, 41, -23, 141, -10, 171, -7,	/* n=64 */
   0, 0, 32, -16, 108, -5, 135, -4};	/* n=128 */
  double conservProb[] =
  {1.00, 0.255, 0.018, 0.008};
  int numConserved[] =
  {64, 138, 30, 12};
  int numDifferent[] =
  {0, 438, 1698, 1716};
  int i, j, k, diff, diffaa, matrixIndex, scoresIndex, changesRead = 0;
  long int changes[4], totalChanges = 0;
  char *codon1, *codon2, label[20];

  fseek (tripfile, -50 * (int) sizeof (long int), SEEK_END);
  while (!feof (tripfile)) {
    fscanf (tripfile, "%s %ld %ld %ld %ld\n", label, &changes[0],
	    &changes[1], &changes[2], &changes[3]);
    if (strcmp (label, "changes:") == 0)
      changesRead = 1;
  }
  if (changesRead == 0) {
    fprintf (stderr, "Changes line not found in triplets file. Aborting\n");
    exit (1);
  }
  else {
    totalChanges = changes[0] + changes[1] + changes[2] + changes[3];
    rewind (tripfile);
  }
  matrix->n = 5;
  matrix->nlabel = (int *) ScoringMatrixCalloc (matrix->n, sizeof (int));
  matrix->score = (int *) ScoringMatrixCalloc (64 * 64 * matrix->n, 
                          sizeof (int));
  matrix->minScore = (int *) ScoringMatrixCalloc (matrix->n, sizeof (int));
  matrix->maxScore = (int *) ScoringMatrixCalloc (matrix->n, sizeof (int));
  matrix->freq = (double *) ScoringMatrixCalloc (64 * 64 * matrix->n, 
                                                 sizeof (double));
  for (i = 0; i < matrix->n; i++) {
    matrix->nlabel[i] = n[i];
    for (j = 0; j < 64; j++) {
      for (k = 0; k < 64; k++) {
	codon1 = codonString (j);
	codon2 = codonString (k);
	diff = codonDiff (codon1, codon2);
	if (translateCodon (codon1, g) != translateCodon (codon2, g))
	  diffaa = 1;
	else
	  diffaa = 0;

	matrixIndex = ((64 * j) + k) + 4096 * i;
	scoresIndex = (i * 8) + (diff * 2) + diffaa;

	matrix->score[matrixIndex] = scores[scoresIndex];
	if (matrix->score[matrixIndex] > matrix->maxScore[i])
	  matrix->maxScore[i] = matrix->score[matrixIndex];
	if (matrix->score[matrixIndex] < matrix->minScore[i])
	  matrix->minScore[i] = matrix->score[matrixIndex];
	if (matrix->score[matrixIndex] >= 0)
	  matrix->freq[matrixIndex] =
	    (conservProb[diff] * changes[diff] / totalChanges) /
	    numConserved[diff];
	else
	  matrix->freq[matrixIndex] = ((1 - conservProb[diff]) * changes[diff] /
				       totalChanges) / numDifferent[diff];
	free (codon1);
	free (codon2);
      }
    }
  }
  return matrix;
}

/* load scoring matrix from file */

ScoringMatrix 
loadMatrix (FILE * file)
{
  ScoringMatrix matrix = (ScoringMatrix) ScoringMatrixCalloc 
                                         (1, sizeof (ScoringMatrixStruct));
  int score, index;
  char dicodon[7];
  double freq;

  matrix->n = 1;
  matrix->nlabel = (int *) ScoringMatrixCalloc (matrix->n, sizeof (int));
  matrix->nlabel[0] = 0;
  matrix->score = (int *) ScoringMatrixCalloc (64 * 64 * matrix->n, 
                                               sizeof (int));
  matrix->minScore = (int *) ScoringMatrixCalloc (matrix->n, sizeof (int));
  matrix->maxScore = (int *) ScoringMatrixCalloc (matrix->n, sizeof (int));
  matrix->freq = (double *) ScoringMatrixCalloc (64 * 64 * matrix->n, 
                                                 sizeof (double));

  while (!feof (file)) {
    fscanf (file, "%s %d %*s %lf\n", dicodon, &score, &freq);
    index = (64 * codonNumber (dicodon)) + codonNumber (&dicodon[3]);
    matrix->score[index] = score;
    matrix->freq[index] = freq;
    if (matrix->score[index] > matrix->maxScore[0])
      matrix->maxScore[0] = matrix->score[index];
    if (matrix->score[index] < matrix->minScore[0])
      matrix->minScore[0] = matrix->score[index];
  }
  fclose (file);
  return matrix;
}

/* return the appropriate comparative score for the two codons */

int
codonScore (ScoringMatrix matrix, int matrixNum, char *codon1,
	    char *codon2)
{
  int cn1 = codonNumber (codon1);
  int cn2 = codonNumber (codon2);

  if ((cn1 >= 0) && (cn2 >= 0))
    return matrix->score[((64 * cn1) + cn2) + 4096 * matrixNum];
  else
    return 0;
}


void
freeScoringMatrix (ScoringMatrix matrix)
{
  free (matrix->nlabel);
  free (matrix->score);
  free (matrix->freq);
  free (matrix);
}

/* wrapper around Calloc to allow error checking */

void *
ScoringMatrixCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for ScoringMatrix\n");
    exit (1);
  }
  else
    return ptr;
}
