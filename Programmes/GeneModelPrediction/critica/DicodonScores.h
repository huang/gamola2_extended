#ifndef __DICODONSCORES_H
#define __DICODONSCORES_H

typedef struct {
  double *score;
  double minScore;
  double maxScore;
  double *freq;
} *DicodonScores, DicodonScoresStruct;

DicodonScores createEmptyDicodonScores ();
DicodonScores loadDicodonScores (FILE *file);
int dicodonNum(char *dicodon);
double scoreDicodon(char *dicodon, DicodonScores dicod);
void freeDicodonScores(DicodonScores dicod);
void *DicodonScoresCalloc (size_t elements, size_t size);

#endif /* __DICODONSCORES_H */
