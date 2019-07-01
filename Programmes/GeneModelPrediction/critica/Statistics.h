#ifndef __STATISTICS_H
#define __STATISTICS_H

typedef struct {
  double *lambda;
  double *k;
  double threshold;
  int strict;
  int sdScan;
  int promScan;
  double alpha;
} *Statistics, StatisticsStruct;

Statistics newStatistics (ScoringMatrix matrix, DicodonScores dicodonScores,
                          int quickStats);
void *StatisticsCalloc (size_t elements, size_t size);
void freeStatistics(Statistics stats);
double KarlinLambdaBis (int low, int high, double *pr, double lambda0);
double KarlinLtoH (double lambda, int low, int high, double *pr);
double KarlinLHtoK (double lambda, double H, int low, int high, double *pr,
		    double score_avg);
double Nlm_Expm1 (register double x);
long Nlm_Gcd (register long a, register long b);

#endif /* __STATISTICS_H */
