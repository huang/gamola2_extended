#ifndef __SEQMAP_H
#define __SEQMAP_H

/* A SeqMap holds a map of a sequence to allow marking of areas that are 
   coding. This is needed for cases in which frequencies of codons coding  
   vs non-coding contexts need to be computed, as well as finding 
   intergenic regions */

typedef struct {
  char *coding;
  int length;
} *SeqMap, SeqMapStruct;

/* creates a new map with the given length */
SeqMap newSeqMap (int length);

/* frees memory used by SeqMap */
void freeSeqMap (SeqMap map);

/* fill a SeqMap from start to end. If AllFrames is 1, mark all frames
   in the region as coding. Otherwise, just mark the translated frame */
void fillSeqMap (SeqMap map, int start, int end, int AllFrames);

int mapPos (SeqMap map, int pos, int dir);
/* wrapper around Calloc */
void *SeqMapCalloc (size_t elements, size_t size);

#endif /* __SEQMAP_H */
