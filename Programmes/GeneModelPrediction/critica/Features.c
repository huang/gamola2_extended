/* 
    Features.c -- Routines for finding SD (RBS) sequences and promoters 
 
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
#include <time.h>
#include "Array.h"
#include "Sequence.h"
#include "GeneticCode.h"
#include "ScoringMatrix.h"
#include "DicodonScores.h"
#include "Statistics.h"
#include "Scores.h"
#include "Regions.h"
#include "Features.h"

void
initCopy (Region * reg, Sequence seq, BonusScores initScores)
{
  strncpy (reg->initSeq, &seq->seq[reg->start], 3);
  reg->initSeq[3] = '\0';
  reg->initScore = findBonus (reg->initSeq, 3, initScores);
}

void
offCopy (Region * reg, BonusScores initScores)
{
  strncpy (reg->initSeq, "off", 3);
  reg->initSeq[3] = '\0';
  reg->initScore = findBonus (reg->initSeq, 3, initScores);
}

/* scans for the best SD sequence for a given initiator at pos */

void
sdScan (Region * reg, Sequence seq, char **matrix,
	BonusScores sdScores)
{
  char sdseq[10] = "RGGRGGTGA";	/* consensus SD sequence */

  int i, j, sdstart = reg->start - 20, match, bestMatch = 0, bestPos = 0;
  int dist = 0, bestDist = 0;

  if (sdstart < 0)
    sdstart = 0;

  for (i = sdstart; i < reg->start - 4; i++) {
    for (j = 0; j < 6; j++) {
      match = findPartialMotif (seq, i + 1, &sdseq[j], matrix);
      dist = reg->start - (i + match) + 1;
      if ((match >= bestMatch) && (dist > 4) && (dist < 14)) {
	bestMatch = match;
	bestPos = i;
	bestDist = dist;
      }
    }
  }
  if (bestMatch > 3) {
    strncpy (reg->sdSeq, &seq->seq[bestPos], bestMatch);
    reg->sdSeq[bestMatch] = '\0';
    reg->sdPos = (short) bestDist;
  }
  reg->sdScore = findBonus (reg->sdSeq, bestMatch, sdScores);
}

/* scans for the best promoter sequence for a given initiator at pos */

#define PROMOTERS 1		/* number of different promoters in list */
#define box35Len 6
#define box10Len 6

void
promoterScan (Region * reg, Sequence seq, char **matrix, 
              BonusScores promScores)
{
  int match35, match10;
  int promoterLen = 0, promStart = 0, len, dist, err10, err35, err;
  int i, j, gap, promDist = 0, promErr=100, searchStart = reg->start - 175;
  char lenStr[3];
  char *box35[PROMOTERS] =
  {"TTGACA"};
  char *box10[PROMOTERS] =
  {"TATAAT"};
  int minGap[PROMOTERS] =
  {16};
  int maxGap[PROMOTERS] =
  {18};

  if (searchStart < 1)
    searchStart = 1;

  for (i = searchStart; i < reg->start; i++) {
    for (j = 0; j < PROMOTERS; j++) {
      match35 = findDegenerateMotif (seq, i + 1, box35[j], matrix, 3, &err35);
      if (match35 > -1) {
	i = match35;
	match10 = findDegenerateMotif (seq, match35 + minGap[j], box10[j],
				       matrix, 3, &err10);
	gap = match10 - (match35 + box35Len - 1);
	len = (match10 + box10Len) - match35;
	dist = reg->start - i - len;
        err = err10 + err35;
	if ((gap >= minGap[j]) && (gap <= maxGap[j]) && (dist > 20) &&
	    (dist < 175) && (err < promErr)) {
	  promoterLen = (match10 + box10Len) - match35 - 1;
	  promStart = i;
	  promDist = dist;
          promErr = err;
	  i += promoterLen;
	}
      }
    }
  }
  if (promoterLen > 0) {
    strncpy (reg->promSeq, &seq->seq[promStart], promoterLen);
    reg->promSeq[promoterLen] = '\0';
    reg->promPos = (short) promDist;
    reg->promErr = (short) promErr;
  }
  sprintf(lenStr,"%d",reg->promErr);
  reg->promScore = findBonus (lenStr, 1, promScores);
}
