#ifndef __SCORES_H
#define __SCORES_H

/* Scores stores the scores calculated from the
   aligned triplets for each sequence position. High scoring regions
   in these scores are considered to be protein coding regions (the
   essence of the CRITICA algorithm). */

typedef struct {
  /* name of contig scores are for */
  char *header;
  /* number of comparative matrices used */
  int n;
  /* array of integer labels for matrices -- 8, 16, etc. */
  int *nlabel;
  /* stores comparative scores for forward strand */
  short **compScores;
  /* stores comparative scores for reverse strand */
  short **antiCompScores;  
  /* length of sequence */
  int seqLength;
} *Scores, ScoresStruct;

/* BonusScores holds scores assigned to the occurrence of a given sequence
   (for example initiator and SD scores) */

typedef struct BonusScoresStruct {
  char *seq;
  double score;
  struct BonusScoresStruct *next;

} *BonusScores, BonusScoresStruct;


/* create a new scores object for a given sequence */

Scores newScores (Sequence seq, ScoringMatrix matrix);

/* compute scores from triplets file and dicodon scores */

Scores loadScores (FILE * file, Sequence seq, Sequence antiSeq,
                   ScoringMatrix matrix);

/* load BonusScores from file */
BonusScores loadBonus (FILE * file);

/* find BonusScore of seq in scores */
double findBonus (char *seq, int len, BonusScores scores);

/* frees memory used by a BonusScores structure */
void freeBonus (BonusScores scores);

void freeScores (Scores compScores);
void *ScoresCalloc (size_t elements, size_t size);

/* prints scores */
void printScores(Scores scores, int dir);

#endif /* __SCORES_H */
