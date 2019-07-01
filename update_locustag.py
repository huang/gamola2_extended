#!/usr/bin/env python


#sed -i 's/"//g' gene_presence_absence.csv
#cp roary/pan_genome_reference.fa gamola/Analysis_Results/HD04_comp1_pan_genome.fa
#cd gamola/Analysis_Results/
#format_fasta_header.py HD04_comp1_pan_genome.fa > HD04_comp1_pan_genome_.fa
#samtools faidx HD04_comp1_pan_genome_.fa
#cp roary/gene_presence_absence.csv gamola/Analysis_Results/


#HD04_01_assem_results/AnnotatedContigs/
#annotated.EFF.reformatted.vcf
#./update_locustag.py HD04_comp1_.fa.gb HD04_comp1_pan_genome_.fa.fai > HD04_comp1__.fa.gb
import sys
import pprint
#from sets import Set
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_dna
from Bio import SeqFeature

pp = pprint.PrettyPrinter(indent=4)
gbk_filename = sys.argv[1]



# import the fai file: mapping the old_id and new_id pan_genome_reference.fa.fai
lines = open(sys.argv[2])
id_names = {}
id_ = 1
for line in lines:
    #newline = line.replace(".", "_")
    tokens = line.split("\t")
    name = tokens[0]
    id_names[id_] = name
    #print("sed -i 's/\/locus_tag=%s/\/locus_tag=%s/g' HD04_comp1_.fa.gb" % (id_, name))
    
    id_ += 1


lines = open(gbk_filename)
for line in lines:
    if "locus_tag=" in line:
        id__ = int(line.split("=")[1].strip().strip("\""))
        print("                     /locus_tag=\"%s\"" % (id_names[id__]))
    else:
        print(line.rstrip())

