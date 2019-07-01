#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Sequence.h"

int
main (int argc, char *argv[])
{
  FILE *seqfile = NULL, *critout = NULL;
  int i, start, end;
  double fraction = 0;
  Sequence seq, antiSeq;
  char *arg, locus[255];
  
  if (argc < 2)
    usage ();

  for (i = 1; i < argc; i++) {
    arg = argv[i];
    if (arg[0] == '-') {
      fprintf (stderr, "%s: option %s invalid\n", argv[0], arg);
      exit (1);
    }
    else if (critout == NULL) {
      critout = fopen (arg, "r");
      if (critout == NULL) {
	fprintf (stderr, "%s: %s not found.\n", argv[0], arg);
	exit (1);
      }
    }
    else if (seqfile == NULL) {
      seqfile = fopen (arg, "r");
      if (seqfile == NULL) {
	fprintf (stderr, "%s: %s not found.\n", argv[0], arg);
	exit (1);
      }
    }
  }
  if ((seqfile == NULL) || (critout == NULL))
    usage ();

  while (!feof (seqfile)) {
    seq = loadNextContig (seqfile);
    antiSeq = complement (seq);
    rewind (critout);
    while (!feof (critout)) {
      fscanf (critout, "%s %d %d%*[^\n]\n", locus, &start, &end);
      if (strcmp (seq->header, locus) == 0) {
	/*fillSeqMap (map, start, end, 0);*/
      }
    }
    freeSequence (seq);
    freeSequence (antiSeq);
  }
  fclose (seqfile);
  fclose (critout);
  free (coding);
  free (noncoding);
  exit (0);
}

void
usage (void)
{
  fprintf (stderr, "Usage: sdscores crit-file seq-file\n");
  exit (1);
}
