#ifndef __SEQUENCE_H
#define __SEQUENCE_H
/* Sequence structures hold a protein or DNA/RNA sequence */

typedef struct {
  int length;			/* length of sequence */
  char *seq;			/* string with sequence data */
  char *header;			/* string for FASTA header */
} *Sequence, SequenceStruct;

/* creates new sequence from string of sequence data and string of header */
Sequence newSequence (char *seqdata, char *header);

/* changes header of sequence structure */
void setHeader (Sequence s, char *newheader);

/* frees memory when done using a Sequence structure */
void freeSequence (Sequence s);

/* prints sequence in FASTA format to file */
void printFasta (Sequence s, FILE * file);

/* returns reverse complement of DNA/RNA sequence */
Sequence complement (Sequence s);

/* returns reverse complement of string -- called by complement */
char *stringComplement (char *seq, int length);

/* return pointer to codon at position pos (first position is 1, as per 
   molecular biology instead of 0 as per C) */
char *codonAt (Sequence s, int pos);

/* returns integer between 0-63 representing codon or -99 if string not a 
   legal codon -- useful for hashing purposes */
int codonNumber (char *codon);

/* returns codon string that would return the number given if used as input
   for codonNumber, above. Remember to free the string when done! */
char *codonString (int num);

/* returns number of positions that are different between two codons
   (0, 1, 2, or 3) */
int codonDiff (char *codon1, char *codon2);

/* returns first position (in molecular biology 1-based numbering) 
   to match to a given motif in a sequence according to the supplied
   pattern matrix.*/
int findMotif(Sequence s, int pos, char *pattern, char **patternMatrix);

/* returns first position (in molecular biology 1-based numbering) 
   to match to a given motif in a sequence according to the supplied
   pattern matrix, with at most err errors */

int  findDegenerateMotif(Sequence s, int pos, char *pattern, 
     char **patternMatrix, int err, int *errfound);

/* returns length of longest match to a given motif at position pos
   in a sequence according to the supplied pattern matrix. */

int findPartialMotif(Sequence s, int pos, char *pattern, char **patternMatrix);

/* returns length of longest match to a given motif at position pos
   in a sequence according to the supplied pattern matrix with at most
   err errors */

int findDegeneratePartialMotif(Sequence s, int pos, char *pattern, 
    char **patternMatrix, int err);
/* returns score of motif at current position from pattern and scores */

int  
scoreMotif(Sequence s, int pos, char *pattern, int *scores,
          char **patternMatrix);

/* returns DNA comparison matrix for use with findMotif, above */
char ** dnaMatrix();

/* returns the most basic matrix for findMotif above -- everything just
   matches itself. */
char ** identityMatrix ();

/* frees memory for matching matrix */
void freeMatchingMatrix(char **matrix);

/* returns G+C content of sequence */
double gcContent(Sequence seq);

/* returns sequence which is a subsequence of the sequence given from 
   start to end (in molecular biology 1-based numbering; in addition, if
   start > end, returns reverse complement of subsequence from end to start 
   -- just like the program "lookat", which not coincidently uses this 
   subroutine */
Sequence lookAt (Sequence s, int start, int end);

/* returns a sequence structure loaded from the next sequence in a FASTA
   file. The function returns NULL if no sequence data was loaded */
Sequence loadNextContig (FILE * file);

/* returns a sequence structure of the sequence with the header equal to
   name. The "seq" field is NULL if no sequence data was loaded */
Sequence loadSpecificContig (char *name, FILE * file);

/* wrappers around Calloc and Realloc to allow error checking */
void *SequenceCalloc (size_t elements, size_t size);
void *SequenceRealloc (void *buffer, size_t size);

#endif /* __SEQUENCE_H */
