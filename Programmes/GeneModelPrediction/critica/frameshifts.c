/* 
    frameshifts.c -- routines to find frameshifts in an array of potential
                     coding regions
 
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
#include "frameshifts.h"

void
frameshiftScan (Array regions, int frameshiftThreshold)
{
  int i, dist, curStart, curEnd;
  int prevStart[3] = {-1,-1,-1};
  int prevEnd[3] = {-1,-1,-1};
  Region *current;

  sortArray (regions, sortRegionsByStart);
  for (i = 0; i < regions->elements; i++) {
    current = (Region *) regions->array[i];
    curStart = realStart (current);
    curEnd = realEnd (current);
    if (prevStart[1+current->dir] > -1) {
      dist = abs (curStart - prevEnd[1+current->dir]);
      if (dist < frameshiftThreshold) {
        fprintf(stderr, 
    "Possible frameshift near %d in %s (would join %s %d %d and %s %d %d)\n",
        prevEnd[1+current->dir], current->name, current->name, 
        prevStart[1+current->dir], prevEnd[1+current->dir], 
        current->name, curStart, curEnd);
      }
    }
    prevStart[1+current->dir] = curStart;
    prevEnd[1+current->dir] = curEnd;
  }
}
