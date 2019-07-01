#ifndef __GENETICCODE_H
#define __GENETICCODE_H

/* Genetic Code structures hold a genetic code table */

typedef struct {
  int num;			/* NCBI genetic code number */
  char *codemap;		/* string containing genetic code */
  short *initiator;            /* array for valid initiator codons*/
} *GeneticCode, GeneticCodeStruct;

/* creates a new GeneticCode structure given an NCBI genetic code number */
GeneticCode newGeneticCode (int num);

/* set the code of a GeneticCode object to a standard NCBI genetic code */
void setCode (GeneticCode g, int num);

/* frees memory used by a GeneticCode structure */
void freeGeneticCode (GeneticCode g);

/* translates the first codon in a string pointed to by *codon */
char translateCodon (char *codon, GeneticCode g);

/* returns 1 if codon pointed to is a valid initiator codon */
int isInitCodon(char *seq, GeneticCode g);

/* returns 1 if codon pointed to is a stop codon */
int isStopCodon(char *seq, GeneticCode g);

/* translates a Sequence structure */
Sequence translate (Sequence s, GeneticCode g);

#endif /* __GENETICCODE_H */
