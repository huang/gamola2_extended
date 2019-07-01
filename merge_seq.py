#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
in_file = open(sys.argv[1],'r')
sequences = SeqIO.parse(in_file, "fasta")
concat=Seq("")
for record in sequences:
    record.seq= record.seq + ("N" * 100)
    concat+=record.seq
print(concat)
