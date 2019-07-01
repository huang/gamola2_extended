#ifndef __REGIONS_H
#define __REGIONS_H

/* Regions structures represent regions believed to be protein-coding */

typedef struct {
  char *name;
  short dir;
  int start;
  int end;
  int seqLength;
  int matrix;
  int n;
  int compScore;
  double dicodonScore;
  double initScore;
  double sdScore;
  short sdPos;
  char sdSeq[11];
  double promScore;
  short promPos;
  short promErr;
  char promSeq[50];
  char initSeq[4];
  double p;
  double k;
  double lambda;
} Region;

Array findRegions (Sequence seq, Sequence antiSeq, Scores scores, 
                   DicodonScores dicod, Statistics stats, GeneticCode g);
void findPeaks(Array array, Sequence seq, short **compScores, int nval, 
               int *nlabel, DicodonScores dicod, Statistics stats, 
               GeneticCode g, int dir);
Region * initRegions(Sequence seq, Statistics stats, int dir, int n, 
                     int label);
int realStart(Region * reg);
int realEnd(Region * reg);
void printRegion (Region *reg, Statistics stats, FILE * file);
void 
extendRegions (Sequence seq, Sequence antiSeq, Array regions, 
               Scores scores, DicodonScores dicod, BonusScores initScores, 
               BonusScores sdScores, BonusScores promScores, Statistics stats,
               GeneticCode g);
Region
recomputeScore (Region reg, short **compScores, Sequence seq,
                DicodonScores dicod, Statistics stats, GeneticCode g);
void
clearScores (Region reg, short **compScores, short **antiCompScores, int n);
void
findInitiator (Region reg, short **compScores, short **antiCompScores,
               int n, Sequence seq, DicodonScores dicod, 
               BonusScores initScores, BonusScores sdScores, 
               BonusScores promScores, Statistics stats, 
               GeneticCode g, char **matrix);
Region
findTerminator (Region reg, Sequence seq, GeneticCode g);
double computeP (Region reg, Statistics stats);
int sortRegionsByP (const void *a, const void *b);
int sortRegionsByStart (const void *a, const void *b);
Region * copyRegion(Region reg);
int regionLen (Region reg);

/* wrappers around Calloc to allow error checking */
void * RegionCalloc (size_t elements, size_t size);

#endif /* __REGIONS_H */
