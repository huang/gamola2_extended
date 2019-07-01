/* 
   Sequence.c -- Data structures for holding a protein or DNA/RNA sequence 

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
#include <ctype.h>
#include <string.h>
#include "Sequence.h"


/* creates new sequence from string of sequence data and string of header */

Sequence
newSequence (char *seqdata, char *header)
{
  Sequence s = (Sequence) SequenceCalloc (1, sizeof (SequenceStruct));

  s->seq = strdup (seqdata);
  s->header = NULL;
  if (s->seq == NULL) {
    fprintf (stderr, "cannot allocate memory for Sequence\n");
    exit (1);
  }
  s->length = strlen (s->seq);
  if (header != NULL)
    setHeader (s, header);
  return s;
}

/* changes header of sequence structure */

void
setHeader (Sequence s, char *newheader)
{
  int i, hlen;

  hlen = (int) strlen(newheader);

  if (s->header != NULL) {
    free (s->header);
  }
  for (i = 0; i < hlen; i++) {
    if (isspace ((int) newheader[i]))
      newheader[i] = '_';
  }
  /* trim excess _'s at end */
  while (newheader[--i]=='_') newheader[i]='\0';
  
  s->header = strdup (newheader);
}

/* frees memory when done using a Sequence structure */

void
freeSequence (Sequence s)
{
  free (s->seq);
  free (s->header);
  free (s);
}

/* prints sequence in FASTA format to file */

void
printFasta (Sequence s, FILE * file)
{
  int i,l;
  fprintf (file, ">%s\n", s->header);
  for (i = 0; i < s->length; i += 60) {
    l = 60; 		/* 60 characters per line */
    if (s->length - i < l) l = s->length - i;
    fprintf (file, "%.*s\n", l, (s->seq)+i);
  }
}


/* returns reverse complement of DNA/RNA sequence */

Sequence
complement (Sequence s)
{
  char *antiseq = stringComplement (s->seq, s->length);
  Sequence s2 = newSequence (antiseq, s->header);
  free (antiseq);
  return s2;
}

/* returns reverse complement of string -- called by complement */

char *
stringComplement (char *seq, int length)
{
  int i;
  char *antiseq;

  antiseq = (char *) SequenceCalloc (length + 1, sizeof (char));
  for (i = 0; i < length; i++) {
    switch (seq[length - i - 1]) {
    case 'A':
      antiseq[i] = 'T';
      break;
    case 'G':
      antiseq[i] = 'C';
      break;
    case 'T':
      antiseq[i] = 'A';
      break;
    case 'C':
      antiseq[i] = 'G';
      break;
    case 'U':
      antiseq[i] = 'A';
      break;
    case 'R':
      antiseq[i] = 'Y';
      break;
    case 'Y':
      antiseq[i] = 'R';
      break;
    case 'W':
      antiseq[i] = 'S';
      break;
    case 'S':
      antiseq[i] = 'W';
      break;
    case 'K':
      antiseq[i] = 'M';
      break;
    case 'M':
      antiseq[i] = 'K';
      break;
    case 'B':
      antiseq[i] = 'V';
      break;
    case 'V':
      antiseq[i] = 'B';
      break;
    case 'D':
      antiseq[i] = 'H';
      break;
    case 'H':
      antiseq[i] = 'D';
      break;
    default:
      antiseq[i] = 'N';
    }
  }
  antiseq[i] = '\0';
  return antiseq;
}

/* return pointer to codon at position pos (first position is 1, as per 
   molecular biology instead of 0 as per C) */

char *
codonAt (Sequence s, int pos)
{
  if ((pos > 0) && (pos < s->length - 1))
    /* only return complete codons, or NULL */
    return &(s->seq[pos - 1]);
  else
    return NULL;
}


/* returns integer between 0-63 representing codon or -99 if string not a 
   legal codon -- useful for hashing purposes */

int
codonNumber (char *codon)
{
  static int fact4[4] =
  {16, 4, 1};
  int i, cnum = 0, value;

  if (codon == NULL)
    return -99;			/* no codon, no number */
  for (i = 0; i < 3; i++) {
    switch (codon[i]) {
    case 'T':
      value = 0;
      break;
    case 'U':
      value = 0;
      break;
    case 'C':
      value = 1;
      break;
    case 'A':
      value = 2;
      break;
    case 'G':
      value = 3;
      break;
    default:
      return -99;
    }
    cnum += fact4[i] * value;
  }
  return cnum;
}


/* returns codon string that would return the number given if used as input
   for codonNumber, above. Remember to free the string when done! */

char *
codonString (int num)
{
  char bases[] = "TCAG";	/* bases assigned to 0, 1, 2, 3 */
  char *codon = (char *) SequenceCalloc (4, sizeof (char));
  char *string = codon;
  int i, rounded, fact4[3] =
  {16, 4, 1};			/* powers of 4 */
  for (i = 0; i < 3; i++) {
    rounded = num / fact4[i];
    *string = bases[rounded];
    num = num - (rounded * fact4[i]);
    string++;
  }
  *string = '\0';
  return codon;
}

/* returns number of positions that are different between two codons
   (0, 1, 2, or 3) */

int
codonDiff (char *codon1, char *codon2)
{
  int differences = 0;
  if (*codon1++ != *codon2++)
    differences++;
  if (*codon1++ != *codon2++)
    differences++;
  if (*codon1 != *codon2)
    differences++;
  return differences;
}


/* returns first position (in molecular biology 1-based numbering) 
   to match to a given motif in a sequence according to the supplied
   pattern matrix. */

int
findMotif (Sequence s, int pos, char *pattern, char **patternMatrix)
{
  int i, j, curr, match, matchPos = -2;
  int succeed = 0, patLen = strlen (pattern);

  for (i = pos - 1; i < s->length - (patLen - 1); i++) {
    succeed = 1;
    for (j = 0; j < patLen; j++) {
      curr = s->seq[i + j];
      match = pattern[j];
      if (patternMatrix[curr][match] == 0) {
	succeed = 0;
	break;
      }
    }
    if (succeed) {
      matchPos = i;
      break;
    }
  }
  return matchPos + 1;
}

/* returns first position (in molecular biology 1-based numbering) 
   to match to a given motif in a sequence according to the supplied
   pattern matrix, with at most err errors */

int
findDegenerateMotif (Sequence s, int pos, char *pattern, char **patternMatrix,
		     int err, int *errfound)
{
  int i, j, curr, match, matchPos = -2;
  int succeed = 0, patLen = strlen (pattern);

  for (i = pos - 1; i < s->length - (patLen - 1); i++) {
    succeed = 1;
    (*errfound) = 0;
    for (j = 0; j < patLen; j++) {
      curr = s->seq[i + j];
      match = pattern[j];
      if (patternMatrix[curr][match] == 0) {
	(*errfound)++;
	if (*errfound > err) {
	  succeed = 0;
	  break;
	}
      }
    }
    if (succeed) {
      matchPos = i;
      break;
    }
  }
  return matchPos + 1;
}

/* returns length of longest match to a given motif at position pos
   in a sequence according to the supplied pattern matrix. */

int
findPartialMotif (Sequence s, int pos, char *pattern, char **patternMatrix)
{
  int i, curr, match;
  int patLen = strlen (pattern);

  for (i = 0; i < patLen; i++) {
    curr = s->seq[pos + i - 1];
    match = pattern[i];
    if (patternMatrix[curr][match] == 0)
      return i;
  }
  return patLen;		/* we didn't return already so we must have matched it all */
}

/* returns length of longest match to a given motif at position pos
   in a sequence according to the supplied pattern matrix with at most
   err errors */

int
findDegeneratePartialMotif (Sequence s, int pos, char *pattern,
			    char **patternMatrix, int err)
{
  int i, curr, match, errors = 0;
  int patLen = strlen (pattern);

  for (i = 0; i < patLen; i++) {
    curr = s->seq[pos + i - 1];
    match = pattern[i];
    if (patternMatrix[curr][match] == 0) {
      errors++;
      if ((errors > err) || (i == 0) || (i == patLen - 1))
	return i;
    }
  }
  return patLen;		/* we didn't return already so we must have matched it all */
}

/* returns score of motif at current position from pattern and scores */

int
scoreMotif (Sequence s, int pos, char *pattern, int *scores,
	    char **patternMatrix)
{
  int i, curr, match, score = 0;
  int patLen = strlen (pattern);

  for (i = 0; i < patLen; i++) {
    if (pos + i >= s->length) {
      score = 0;
      break;
    }
    curr = s->seq[pos + i - 1];
    match = pattern[i];
    if (patternMatrix[curr][match] > 0) {
      score += scores[i];
    }
  }
  return score;
}


/* returns DNA comparison matrix for use with findMotif, above */

char **
dnaMatrix ()
{
  char **matrix = identityMatrix ();
  int i;

  /* R = A|G */
  matrix['R']['A'] = 1;
  matrix['A']['R'] = 1;
  matrix['R']['G'] = 1;
  matrix['G']['R'] = 1;

  /* Y = C|T */
  matrix['Y']['C'] = 1;
  matrix['C']['Y'] = 1;
  matrix['Y']['T'] = 1;
  matrix['T']['Y'] = 1;

  /* W = A|T */
  matrix['W']['A'] = 1;
  matrix['A']['W'] = 1;
  matrix['W']['T'] = 1;
  matrix['T']['W'] = 1;

  /* S = G|C */
  matrix['S']['G'] = 1;
  matrix['G']['S'] = 1;
  matrix['S']['C'] = 1;
  matrix['C']['S'] = 1;

  /* N matches anything */
  for (i = 0; i < 255; i++) {
    matrix['N'][i] = 1;
    matrix[i]['N'] = 1;
  }

  return matrix;
}

/* returns the most basic matrix for findMotif above -- everything just
   matches itself. */

char **
identityMatrix ()
{
  int i;
  char **matrix;

  matrix = (char **) SequenceCalloc (255, sizeof (char *));

  for (i = 0; i < 255; i++)
    matrix[i] = (char *) SequenceCalloc (255, sizeof (char));
  for (i = 0; i < 255; i++)
    matrix[i][i] = 1;
  return matrix;
}

/* frees memory for matching matrix */

void
freeMatchingMatrix (char **matrix)
{
  int i;

  for (i = 0; i < 255; i++)
    free (matrix[i]);
  free (matrix);
}

/* returns G+C content of sequence */

double
gcContent (Sequence seq)
{
  int i;
  double val = 0;

  for (i = 0; i < seq->length; i++)
    if ((seq->seq[i] == 'G') || (seq->seq[i] == 'C'))
      val++;
  val /= seq->length;

  return val;
}

/* returns sequence which is a subsequence of the sequence given from 
   start to end (in molecular biology 1-based numbering; in addition, if
   start > end, returns reverse complement of subsequence from end to start 
   -- just like the program "lookat", which not coincidently uses this 
   subroutine */

Sequence
lookAt (Sequence s, int start, int end)
{
  int sublen;
  char *str, hstring[550];
  Sequence subseq, antisubseq;

  /* bounds checking */
  if (start > s->length)
    start = s->length;
  if (start < 1)
    start = 1;
  if (end > s->length)
    end = s->length;
  if (end < 1)
    end = 1;

  /* header of subsequence includes position information */
  sprintf (hstring, "%s_%d_%d", s->header, start, end);

  if (start <= end) {		/* start is less than end -- normal case */
    sublen = end + 1 - start;
    str = (char *) SequenceCalloc (sublen + 1, sizeof (char));
    strncpy (str, &(s->seq[start - 1]), sublen);
    subseq = newSequence (str, hstring);
  }
  else {	  /* start is greater than end -- need reverse complement */
    sublen = start + 1 - end;
    str = (char *) SequenceCalloc (sublen + 1, sizeof (char));
    strncpy (str, &(s->seq[end - 1]), sublen);
    subseq = newSequence (str, hstring);
    antisubseq = complement (subseq);
    freeSequence (subseq);
    subseq = antisubseq;
  }
  free (str);
  return subseq;
}


#define BUFSIZE 10240		/* constant amount of reallocation used in 
				   loadNextContig */

/* returns a sequence structure loaded from the next sequence in a FASTA
   file. The "seq" field is NULL if no sequence data was loaded */

Sequence
loadNextContig (FILE * file)
{
  int i = 0, bufsize = BUFSIZE, newchar;
  char header[500] = "\0";
  char *buffer = (char *) SequenceCalloc (BUFSIZE, sizeof (char));
  Sequence s = (Sequence) SequenceCalloc (1, sizeof (SequenceStruct));

  s->seq = NULL;
  newchar = getc (file);
  if (newchar != '>')
    ungetc (newchar, file);	/* oops -- this isn't a FASTA header! */
  else
    fscanf (file, "%499[^\n]\n", header);

  while (((newchar = getc (file)) != EOF) && (newchar != '>')) {
    if (i + 1 > bufsize) {	/* If buffer is full, resize it. */
      bufsize += BUFSIZE;
      buffer = (char *) SequenceRealloc (buffer, bufsize);
    }
    if ((isalpha (newchar)) || (newchar == '*')) {  /* only add valid chars */
      buffer[i] = (char) toupper (newchar);
      i++;
    }
  }
  buffer[i] = '\0';		/* Tack on a null-terminator. */
  if (newchar == '>')
    ungetc (newchar, file);	/* put back '>' of next seq */
  if (newchar == EOF && *buffer == '\0') {
    s = NULL;
  }
  else {
    s = newSequence (buffer, header);
  }
  free (buffer);
  return s;
}


/* returns a sequence structure of the sequence with the header equal to
   name. The function returns  NULL if no sequence data was loaded */

Sequence
loadSpecificContig (char *name, FILE * file)
{
  Sequence s = NULL;
  int i;
  
  for (i = 0; i < (int) strlen (name); i++) {
    if (isspace ((int) name[i]))
      name[i] = '_';
  }
   
  while (!feof (file)) {
    s = loadNextContig (file);
    if (s != NULL) {
      if (strcmp (name, s->header) == 0)
	break;
      else {
	freeSequence (s);
	s = NULL;
      }
    }
  }
  return s;
}

/* wrappers around Calloc and Realloc to allow error checking */

void *
SequenceCalloc (size_t elements, size_t size)
{
  void *ptr;
  ptr = calloc (elements, size);
  if (ptr == NULL) {
    fprintf (stderr, "cannot allocate memory for Sequence\n");
    exit (1);
  }
  else
    return ptr;
}

void *
SequenceRealloc (void *buffer, size_t size)
{
  buffer = realloc (buffer, size);
  if (buffer == NULL) {
    fprintf (stderr, "cannot reallocate memory for Sequence\n");
    exit (1);
  }
  else
    return buffer;
}
