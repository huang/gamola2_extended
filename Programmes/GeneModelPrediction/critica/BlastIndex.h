#ifndef __BLASTINDEX_H
#define __BLASTINDEX_H

typedef struct BlastIndexNode {
  char *name;
  double p;
  int start;
  int end;
  long matchPos;
  struct BlastIndexNode *next;
} BlastIndexNode, *BlastIndex;

BlastIndex newBlastIndex (FILE * file);
BlastIndex nextBlastNode (BlastIndex index);
void freeBlastIndex (BlastIndex index);
void *BlastIndexCalloc (size_t elements, size_t size);

#endif /* __BLASTINDEX_H */
