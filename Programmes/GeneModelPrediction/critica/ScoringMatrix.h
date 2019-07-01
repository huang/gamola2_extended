#ifndef __SCORINGMATRIX_H
#define __SCORINGMATRIX_H

/* a ScoringMatrix structure holds a lookup table of comparative scores 
   for every pair of codons which CompScores uses to calculate the comparative
   score for each position */

typedef struct {
  int *score;
  int *minScore;
  int *maxScore;
  double *freq;
  int *nlabel;
  int n;
} *ScoringMatrix, ScoringMatrixStruct;

/* creates the default set of scoring matrices as described in the CRITICA
   paper */
ScoringMatrix createDefaultScoringMatrix (GeneticCode g, FILE * tripfile);

/* load scoring matrix from file */

ScoringMatrix loadMatrix(FILE *file);

/* return the appropriate comparative score for the two codons */
int codonScore (ScoringMatrix matrix, int matrixNum, char *codon1,
		char *codon2);

void freeScoringMatrix (ScoringMatrix matrix);

void *ScoringMatrixCalloc (size_t elements, size_t size);

#endif /* __SCORINGMATRIX_H */
