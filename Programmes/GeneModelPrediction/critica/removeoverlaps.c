/* 
    removeoverlaps.c -- Program for removing coding regions inconsistent
                        with coding regions with more support
 
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
#include "Array.h"
#include "Sequence.h"
#include "SeqMap.h"

typedef struct Call {
  char *locus;
  char *line;
  int start;
  int end;
  double p;
} Call;

void usage (void);
int computeOverlap(Call *call, SeqMap map);
int compareCall (const void *a, const void *b);
void freeCall (void *a);


int
main (int argc, char *argv[])
{
  FILE *seqfile, *cdsfile;
  int i, percent, start, end, overlap;
  double p;
  char locus[500],line[500];
  Sequence seq, antiSeq;
  SeqMap map;
  Call *call;
  Array calls=newArray(5000);

  if (argc != 4)
    usage ();

  seqfile=fopen(argv[1],"r");
  if (seqfile == NULL) {
    fprintf(stderr,"%s: can't open %s\n",argv[0],argv[1]);
    exit(1);
  }
  cdsfile=fopen(argv[2],"r");
  if (cdsfile == NULL) {
    fprintf(stderr,"%s: can't open %s\n",argv[0],argv[2]);
    exit(1);
  }
  percent=atoi(argv[3]);
  
  while (!feof (cdsfile)) {
    line[0]='\0';
    fscanf(cdsfile, "%499[^\n]\n", line);
    sscanf (line, "%s %d %d %le %*[^\n]", locus, &start, &end, &p);
    if (start==0) continue;
    call= (Call *) SeqMapCalloc(1, sizeof(Call));
    call->locus=strdup(locus);
    call->line=strdup(line);
    call->start=start;
    call->end=end;
    call->p=p;
    addToArray(calls, call);
  }
  sortArray(calls, compareCall);
  while (!feof (seqfile)) {
    seq = loadNextContig (seqfile);
    antiSeq = complement (seq);
    map = newSeqMap (seq->length);
    for(i=0;i<calls->elements;i++) {
       call = (Call *) calls->array[i];
       if (strcmp (seq->header, call->locus) == 0) {
         overlap=computeOverlap(call,map);
         if (overlap>=percent) 
           fprintf(stderr, "%s\n", call->line);
	 else
           fprintf(stdout, "%s\n", call->line);
         fillSeqMap (map, call->start, call->end, 1);
       }
     }
    freeSequence(seq);
    freeSequence(antiSeq);
    freeSeqMap(map);
  }
  freeArray(calls, freeCall);
  fclose (cdsfile); 
  fclose (seqfile);
  exit (0);
}

int computeOverlap(Call *call, SeqMap map) 
{
  int overlap=0;
  int i, len=1+abs(call->start-call->end);

  if (call->start < call->end) {
    for(i=call->start;i<=call->end;i++) {
      if (mapPos(map,i-1,1)) overlap++;
    }
  }
  else {
    for(i=call->start-1;i>call->end;i--) {
      if (mapPos(map,i-1,1)) overlap++;
    }
  }
  return (overlap*100)/len;
}

int
compareCall (const void *a, const void *b)
{
  const Call *call1 = (Call *) a, *call2 = (Call *) b;
  if (call1->p > call2->p)
    return 1;
  else if (call1->p < call2->p)
    return -1;
  else
    return 0;
}

void freeCall (void *a) 
{
  Call *call = (Call *) a;
  
  free(call->locus);
  free(call->line);
  free(call);
}

void
usage (void)
{
  fprintf (stderr, 
    "Usage: removeoverlaps seq-file cds-file percent-allowed\n");
  exit (1);
}
