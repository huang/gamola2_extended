/* 
    Statistics.c -- Routines for computing the significance of a coding
                    region
 
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
#include <math.h>
#include "Sequence.h"
#include "GeneticCode.h"
#include "ScoringMatrix.h"
#include "DicodonScores.h"
#include "Statistics.h"


Statistics
newStatistics (ScoringMatrix matrix, DicodonScores dicodonScores,
	       int quickStats)
{
  Statistics stats = (Statistics) StatisticsCalloc (1, 
                                  sizeof (StatisticsStruct));
  int i, j, k, score;
  double roughLambda = 0.015;
  int min, max;
  double *pr, h, prob, totalProb, averageScore;

  stats->threshold = 1e-4;
  stats->strict = 0;
  stats->alpha = 0.8;
  stats->sdScan = 1;
  stats->promScan = 0;
  stats->lambda = (double *) StatisticsCalloc (matrix->n, sizeof (double));
  stats->k = (double *) StatisticsCalloc (matrix->n, sizeof (double));
  for (i = 0; i < matrix->n; i++) {
    if (quickStats) {
      stats->lambda[i] = 0.015;
      stats->k[i] = 0.2;
    }
    else {
      min = matrix->minScore[i] + (int) (dicodonScores->minScore / 
                                         roughLambda);
      max = matrix->maxScore[i] + (int) (dicodonScores->maxScore / 
                                         roughLambda);
      pr = (double *) StatisticsCalloc (1 + max - min, sizeof (double));
      totalProb = 0;
      averageScore = 0;
      for (j = 0; j < 4096; j++) {
	for (k = 0; k < 4096; k++) {
	  score = matrix->score[(4096 * i) + j] +
	    (int) (dicodonScores->score[k] / roughLambda);
	  prob = matrix->freq[(4096 * i) + j] * dicodonScores->freq[k];
	  pr[score - min] += prob;
	  totalProb += prob;
	  averageScore += (score * prob);
	}
      }
      if ((totalProb < .99) || (totalProb > 1.01)) {
	fprintf (stderr, "Total probability %e should be 1\n", totalProb);
	exit (1);
      }
      if (averageScore > 0) {
	fprintf (stderr, 
          "Help! Positive Average Score %8.3e!\n", averageScore);
	exit (1);
      }
      stats->lambda[i] = KarlinLambdaBis (min, max, pr, 0.001);
      h = KarlinLtoH (stats->lambda[i], min, max, pr);
      stats->k[i] = KarlinLHtoK (stats->lambda[i], h, min, max, pr, 
                                 averageScore);
      free (pr);
    }
  }
  return stats;
}

void 
freeStatistics (Statistics stats)
{
  free (stats->lambda);
  free (stats->k);
  free (stats);
}

/* wrapper around Calloc to allow error checking */

void *
StatisticsCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for Statistics\n");
    exit (1);
  }
  else
    return ptr;
}


/* Karlin-Altschul Max Segment statistics routines taken from BLAST, and were
   presumably written by Steven Altschul */

/*
 * KarlinLambdaBis
 * 
 * Calculate Lambda using the bisection method (slow).
 */

double
KarlinLambdaBis (int low, int high, double *pr, double lambda0)
{
  double lambda, up, newval;
  int i, j;
  register double sum, x0, x1;


  up = lambda0;
  for (lambda = 0.;;) {
    up *= 2;
    x0 = exp ((double) up);
    x1 = pow ((double) x0, low - 1);
    if (x1 > 0.) {
      for (sum = 0., i = low; i <= high; ++i)
	sum += pr[i - low] * (x1 *= x0);
    }
    else {
      for (sum = 0., i = low; i <= high; ++i)
	sum += pr[i - low] * exp (up * i);
    }
    if (sum >= 1.0)
      break;
    lambda = up;
  }

  for (j = 0; j < 20; ++j) {
    newval = (lambda + up) / 2.;
    x0 = exp ((double) newval);
    x1 = pow ((double) x0, low - 1);
    if (x1 > 0.) {
      for (sum = 0., i = low; i <= high; ++i)
	sum += pr[i - low] * (x1 *= x0);
    }
    else {
      for (sum = 0., i = low; i <= high; ++i)
	sum += pr[i - low] * exp (newval * i);
    }
    if (sum > 1.0)
      up = newval;
    else
      lambda = newval;
  }
  return (lambda + up) / 2.;
}


/*
 * KarlinLtoH
 * 
 * Calculate H, the relative entropy of the p's and q's
 */
double
KarlinLtoH (double lambda, int low, int high, double *pr)
{
  register int i;
  register double av, etolam, etolami;


  etolam = exp ((double) lambda);
  etolami = pow ((double) etolam, low - 1);
  if (etolami > 0.) {
    for (av = 0., i = low; i <= high; ++i)
      av += pr[i - low] * i * (etolami *= etolam);
  }
  else {
    for (av = 0., i = low; i <= high; ++i)
      av += pr[i - low] * i * exp (lambda * i);
  }
  return lambda * av;
}

#define DIMOFP0(iter,range)	(iter*range + 1)

double
KarlinLHtoK (double lambda, double H, int low, int high, double *pr, double score_avg)
{

  double *P0 = NULL;

  double K;			/* local copy of K */
  double ratio;
  int i, j;
  int range, lo, hi, first, last;
  register double sum;
  double Sum, av, oldsum, oldsum2;
  int iter;
  double sumlimit, x;
  double *p, *ptrP, *ptr1, *ptr2, *ptr1e;
  double etolami, etolam;

  if (low == -1 || high == 1) {
    if (low == -1 && high == 1) {
      x = (pr[low - low] - pr[high - low]);
      K = x * x / pr[low - low];
      return K;
    }
    av = H / lambda;
    if (high == 1)
      K = av;
    else
      K = (score_avg * score_avg) / av;
    K *= -Nlm_Expm1 (-lambda);
    return K;
  }
  if (high == -low) {
    /* see if there is only one positive score observed, high */
    for (first = 1; first < high; ++first) {
      if (pr[first - low] > 0.)
	break;
    }
    if (first == high) {
      av = H / (lambda *= high);
      K = av * -Nlm_Expm1 (-lambda);
      return K;
    }
  }
  sumlimit = 0.01;
  range = high - low;
  iter = 17;

  P0 = (double *) calloc (DIMOFP0 (iter, range), sizeof (double));
  if (P0 == NULL)
    return -1.;
  av = H / lambda;
  etolam = exp ((double) lambda);
  Sum = 0.;
  lo = hi = 0;
  p = &pr[low - low];
  P0[0] = sum = oldsum = oldsum2 = 1.;
  for (j = 0; j < iter && sum > sumlimit; Sum += sum /= ++j) {
    first = last = range;
    lo += low;
    hi += high;
    for (ptrP = P0 + (hi - lo); ptrP >= P0; *ptrP-- = sum) {
      ptr1 = ptrP - first;
      ptr1e = ptrP - last;
      ptr2 = p + first;
      for (sum = 0.; ptr1 >= ptr1e;)
	sum += *ptr1-- * *ptr2++;
      if (first)
	--first;
      if (ptrP - P0 <= range)
	--last;
    }
    etolami = pow ((double) etolam, lo - 1);
    for (sum = 0., i = lo; i != 0; ++i) {
      etolami *= etolam;
      sum += *++ptrP * etolami;
    }
    for (; i <= hi; ++i)
      sum += *++ptrP;
    oldsum2 = oldsum;
    oldsum = sum;
  }

  /* Terms of geometric progression added for correction */
  ratio = oldsum / oldsum2;
  if (ratio >= (1.0 - sumlimit * 0.001)) {
    K = -1.;
    goto CleanUp;
  }
  sumlimit *= 0.01;
  while (sum > sumlimit) {
    oldsum *= ratio;
    Sum += sum = oldsum / ++j;
  }

  for (i = 1, j = -low; i <= range && j > 1; ++i)
    if (p[i])
      j = Nlm_Gcd (j, i);

  if (j * etolam > 0.05) {
    etolami = pow ((double) etolam, -j);
    K = j * exp ((double) -2.0 * Sum) / (av * (1.0 - etolami));
  }
  else
    K = -j * exp ((double) -2.0 * Sum) / (av * Nlm_Expm1 (-j * (double) lambda));

CleanUp:
  if (P0 != NULL)
    free (P0);
  return K;
}

double
Nlm_Expm1 (register double x)
{
  register double absx;

  if ((absx = fabs (x)) > .33)
    return exp (x) - 1.;

  if (absx < 1.e-16)
    return x;

  return x * (1. + x * (0.5 + x * (1. / 6. + x * (1. / 24. + x * (1. /
			 120. + x * (1. / 720. + x * (1. / 5040. + x * (1. /
		       40320. + x * (1. / 362880. + x * (1. / 3628800. + x *
				(1. / 39916800. + x * (1. / 479001600. + x /
						    6227020800.))))))))))));
}

long
Nlm_Gcd (register long a, register long b)
{
  register long c;

  b = abs (b);
  if (b > a)
    c = a, a = b, b = c;

  while (b != 0) {
    c = a % b;
    a = b;
    b = c;
  }
  return a;
}
