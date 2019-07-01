#ifndef __CODONNODE_H
#define __CODONNODE_H

typedef struct CodonNode {
  unsigned char codon;
  short num;
  struct CodonNode *next;
} CodonNode;

typedef struct CodonListStruct {
  char *name;
  int length;
  CodonNode **list;
} *CodonList, CodonListStruct;


CodonList newCodonList (char *name, int length);
void addCodon (CodonList list, char *codon, int pos);
void printCodonList (CodonList list, FILE * file);
void freeCodonList (CodonList list);
void *CodonListCalloc (size_t elements, size_t size);

#endif /* __CODONNODE_H */
