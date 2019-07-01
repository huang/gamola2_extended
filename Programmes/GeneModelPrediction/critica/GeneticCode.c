/* 
    GeneticCode.c -- Data structure for holding a genetic code table 
 
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


#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Sequence.h"
#include "GeneticCode.h"

/* creates a new GeneticCode structure given an NCBI genetic code number */

GeneticCode 
newGeneticCode (int num)
{
  GeneticCode g= (GeneticCode) SequenceCalloc(1, sizeof(GeneticCodeStruct));
  g->codemap = NULL;
  setCode (g, num);
  return g;
}

/* frees memory used by a GeneticCode structure */

void 
freeGeneticCode (GeneticCode g)
{
  free (g->codemap);
  free (g->initiator);
  free(g);
}

/* translates the first codon in a string pointed to by *codon */

char 
translateCodon (char *codon, GeneticCode g)
{
  int cn = codonNumber (codon);
  if (cn < 0)
    return 'X';
  else
    return g->codemap[cn];
}

/* returns 1 if codon pointed to is a valid initiator codon */

int isInitCodon(char *seq, GeneticCode g) 
{
  int cn = codonNumber (seq);
  if (cn < 0) 
    return 0;
  else 
   return g->initiator[cn];
}

/* returns 1 if codon pointed to is a stop codon */

int isStopCodon(char *seq, GeneticCode g) 
{
  int cn = codonNumber (seq);
  if (cn < 0) 
    return 0;
  else 
    if (g->codemap[cn]=='*')
      return 1;
  else
    return 0; 
}


/* translates a Sequence structure */

Sequence 
translate (Sequence s, GeneticCode g)
{
  int i;
  char *protstr = (char *) SequenceCalloc ((s->length / 3) + 1, sizeof (char));
  char *p = protstr, hstring[550];
  Sequence protseq;

  for (i = 0; i < s->length - 2; i += 3) {
    *p++ = translateCodon (&s->seq[i], g);
  }
  *p = '\0';
  sprintf (hstring, "%s (translated)", s->header);
  protseq = newSequence (protstr, hstring);
  free (protstr);
  return protseq;
}

/* set the code of a GeneticCode object to a standard NCBI genetic code */

void 
setCode (GeneticCode g, int num)
{
  if (g->codemap != NULL) {
    free(g->codemap);
    free(g->initiator);
  }
  g->num=num;

  switch (num) {
  case 1:			/* Standard */
    g->codemap = strdup(
       "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 2:			/* Vertebrate Mitochondrial */
    g->codemap = strdup(
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 3:			/* Yeast Mitochondrial */
    g->codemap = strdup(
	"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 4:			/* Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma */
    g->codemap = strdup(
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 5:			/* Invertebrate Mitochondrial */
    g->codemap = strdup(
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 6:			/* Ciliate Macronuclear and Dasycladacean */
    g->codemap = strdup(
	"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 9:			/* Echinoderm Mitochondrial */
    g->codemap = strdup(
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 10:			/* Euplotid Nuclear */
    g->codemap = strdup(
	"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 11:			/* Bacterial */
    g->codemap = strdup(
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    g->initiator[codonNumber("TTG")]=1;
    g->initiator[codonNumber("GTG")]=1;
    break;
  case 12:			/* Alternative Yeast Nuclear */
    g->codemap = strdup(
	"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 13:			/*  Ascidian Mitochondrial */
    g->codemap = strdup(
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 14:			/* Flatworm Mitochondrial */
    g->codemap = strdup(
	"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  case 15:			/* Blepharisma Macronuclear */
    g->codemap = strdup(
	"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    g->initiator= (short *) SequenceCalloc(64, sizeof(short));
    g->initiator[codonNumber("ATG")]=1;
    break;
  default:
    fprintf (stderr, "No such genetic code %d\n", num);
    exit (1);
  }
}
