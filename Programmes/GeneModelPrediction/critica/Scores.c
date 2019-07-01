/* 
    Scores.c -- Data structure for holding comparative scores
 
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
#include <string.h>
#include "Sequence.h"
#include "GeneticCode.h"
#include "ScoringMatrix.h"
#include "DicodonScores.h"
#include "Scores.h"


/* create a new scores object for a given sequence */

Scores
newScores (Sequence seq, ScoringMatrix matrix)
{
  int i;
  Scores scores= (Scores) ScoresCalloc(1, sizeof(ScoresStruct));

  scores->seqLength = seq->length;
  scores->header = strdup (seq->header);
  scores->compScores = (short **) ScoresCalloc (matrix->n, sizeof (short *));
  scores->antiCompScores = (short **) ScoresCalloc (matrix->n, 
                                      sizeof (short *));
  scores->nlabel = (int *) ScoresCalloc (matrix->n, sizeof (int));
  scores->n = matrix->n;
  for (i = 0; i < matrix->n; i++) {
    scores->compScores[i] = (short *) ScoresCalloc (seq->length, 
                                      sizeof (short));
    scores->antiCompScores[i] = (short *) ScoresCalloc (seq->length, 
                                          sizeof (short));
    scores->nlabel[i] = matrix->nlabel[i];
  }
  return scores;
}


/* compute scores from triplets file and dicodon scores */
Scores
loadScores (FILE * file, Sequence seq, Sequence antiSeq, ScoringMatrix matrix)
{
  char line[1000], orig[1000],*tok, *codon1, codon2[4]; 
  char *antiCodon1, *antiCodon2;
  short *sites;
  Scores scores;
  int i, j, num, antiNum, count, score, antiScore;
  long pos;

  scores = newScores (seq, matrix);
  sites = (short *) ScoresCalloc (seq->length, sizeof (short));

  fscanf(file,"%s\n",line);
  if (strcmp(seq->header,line)!=0) {
    fprintf(stderr,"critica: fatal error -- was looking for %s in triplets but found %s instead\n",seq->header,line);
   exit(1);   
  }
  while (!feof (file)) {
    pos=ftell(file);  /* remember position of start of line */
    fscanf (file, "%999[^\n]\n", line);
    strcpy(orig, line);
    tok = strtok (line, " ");
    if ((strcmp(tok, orig)==0) || (strcmp(tok,"changes:")==0)){
      /* no more triplets for sequence -- we're done */
      fseek (file, pos, SEEK_SET);  
      
      break;
    }
    num = atoi(tok);
    antiNum = seq->length - num - 1;
    codon1 = codonAt (seq, num);
    antiCodon1 = codonAt (antiSeq, antiNum);
    tok = strtok (NULL, " ");

    while (tok != NULL) {
      strcpy (codon2, tok);
      antiCodon2 = stringComplement (codon2, 3);
      count = atoi (strtok (NULL, " "));
      for (i = 0; i < matrix->n; i++) {
	score = codonScore (matrix, i, codon1, codon2);
	scores->compScores[i][num - 1] += count * score;
	antiScore = codonScore (matrix, i, antiCodon1, antiCodon2);
	scores->antiCompScores[i][antiNum - 1] += (short) (count * antiScore);
      }
      sites[num - 1] += count;
      free (antiCodon2);
      tok = strtok (NULL, " ");
    }
  }
  for (j = 0; j < scores->seqLength - 2; j++) {
    if (sites[j] > 1) {
      antiNum = seq->length - (j+1) - 1;
      for (i = 0; i < matrix->n; i++) {
	scores->compScores[i][j] /= sites[j];
	scores->antiCompScores[i][antiNum - 1] /= sites[j];
      }
    }
  }
  free (sites);
  return scores;
}

/* load BonusScores from file */

BonusScores
loadBonus (FILE * file)
{
  BonusScores start = NULL, current, previous = NULL;
  char seq[255];
  double score;

  while (!feof (file)) {
    fscanf (file, "%s %*d %*d %*f %*f %lf\n", seq, &score);
    current = (BonusScores) ScoresCalloc (1, sizeof (BonusScoresStruct));
    current->seq = strdup (seq);
    current->score = score;
    current->next = NULL;
    if (previous != NULL)
      previous->next = current;
    else
      start = current;
    previous = current;
  }
  return start;
}

/* find BonusScore of seq in scores */
double
findBonus (char *seq, int len, BonusScores scores)
{
  BonusScores current = scores;

  while (current != NULL) {
    if (strncmp (seq, current->seq, len) == 0)
      return current->score;
    current = current->next;
  }
  return 0;
}

/* frees memory used by a BonusScores structure */
void
freeBonus (BonusScores scores)
{
  BonusScores current = scores, next;
  while (current != NULL) {
    next = current->next;
    free (current->seq);
    free (current);
    current = next;
  }
}


/* frees memory used by a Scores structure */
void
freeScores (Scores scores)
{
  int i;
  
  for (i = 0; i < scores->n; i++) {
    free (scores->compScores[i]);
    free (scores->antiCompScores[i]);
  }
  free (scores->compScores);
  free (scores->antiCompScores);
  free (scores->nlabel);
  free (scores->header);
  free (scores);
}

/* prints scores */

void 
printScores (Scores scores, int dir)
{
  int i;

  for (i = 0; i < scores->seqLength; i++) {
    if (dir == 1)
      fprintf (stderr, "%5d %5d\n", i + 1, scores->compScores[0][i]);
    else
      fprintf (stderr, "%5d %5d\n", i + 1, scores->antiCompScores[0][i]);
  }
}

/* wrapper around Calloc to allow error checking */

void *
ScoresCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for Scores\n");
    exit (1);
  }
  else
    return ptr;
}
