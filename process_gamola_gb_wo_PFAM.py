#!/usr/bin/env python


#ftp://ftp.ncbi.nlm.nih.gov/blast/db/

#./extract_product_from_gbk.py results/featureCounts_foldchange3/adeIJ_ko_vs_wt_background.txt BIT33__.fa.gb > adeIJ_ko_vs_wt_annotated_degenes.csv
#./extract_product_from_gbk.py results/featureCounts_foldchange3/craA_ko_vs_wt_background.txt BIT33__.fa.gb > craA_ko_vs_wt_annotated_degenes.csv
#./extract_product_from_gbk.py results/featureCounts_foldchange3/adeIJ_ko_vs_craA_ko_background.txt BIT33__.fa.gb > adeIJ_ko_vs_craA_ko_annotated_degenes.csv
#"BIT33_RS06380"


#~/Tools/csv2xls-0.4/csv_to_xls.py craA_ko_vs_wt.csv adeIJ_ko_vs_wt.csv adeIJ_ko_vs_craA_ko.csv -d$',' -o annotated_degenes.xls
#~/Tools/csv2xls-0.4/csv_to_xls.py craA_ko_vs_wt.csv -d$',' -o annotated_degenes.xls

#sed -i 's/"//g' gene_presence_absence.csv
#cp roary/pan_genome_reference.fa gamola/Analysis_Results/HD04_comp1_pan_genome.fa
#cd gamola/Analysis_Results/
#format_fasta_header.py HD04_comp1_pan_genome.fa > HD04_comp1_pan_genome_.fa
#samtools faidx HD04_comp1_pan_genome_.fa
#cp roary/gene_presence_absence.csv gamola/Analysis_Results/



## REFORMAT the fa.gb file
#cp GAMOLA2/Consolidated_results/BIT33.fa/BIT33.fa.gb ./
#awk '{print $0 "ORIGIN"> "file" NR}' RS='ORIGIN'  BIT33.fa.gb
#sed -i 's/</ /g' file2
#sed -i 's/^ //g' file2
#cat file1 file2 > BIT33.fa.gb_.fa.gb
#samtools faidx results/reference_genome/GCF_002811175.fna
#update_locustag.py example_.fa.gb results/reference_genome/GCF_002811175.fna.fai > BIT33.fa.gb__.fa.gb
## manually clean the empty lines after 'ORIGIN' and the last line containing 'ORIGIN'



#HD04_01_assem_results/AnnotatedContigs/
#annotated.EFF.reformatted.vcf
#./extract_product_from_gbk.py adeIJ_ko_vs_wt_background.txt BIT33__.fa.gb > annotated_degenes.csv
#delete the empty three columns.
import sys
import pprint
#from sets import Set
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_dna
from Bio import SeqFeature

pp = pprint.PrettyPrinter(indent=4)

input_for_reformatting = sys.argv[1]
gbk_filename = sys.argv[2]
#ncbi_gbk_filename = sys.argv[3]
#tigrfam_gbk_filename = sys.argv[4]



'''
# Code  Name
J       Translation, ribosomal structure and biogenesis
A       RNA processing and modification
K       Transcription
L       Replication, recombination and repair
B       Chromatin structure and dynamics
D       Cell cycle control, cell division, chromosome partitioning
Y       Nuclear structure
V       Defense mechanisms
T       Signal transduction mechanisms
M       Cell wall/membrane/envelope biogenesis
N       Cell motility
Z       Cytoskeleton
W       Extracellular structures
U       Intracellular trafficking, secretion, and vesicular transport
O       Posttranslational modification, protein turnover, chaperones
X       Mobilome: prophages, transposons
C       Energy production and conversion
G       Carbohydrate transport and metabolism
E       Amino acid transport and metabolism
F       Nucleotide transport and metabolism
H       Coenzyme transport and metabolism
I       Lipid transport and metabolism
P       Inorganic ion transport and metabolism
Q       Secondary metabolites biosynthesis, transport and catabolism
R       General function prediction only
S       Function unknown
'''

locustag_gene = {}
locustag_product = {}
locustag_note = {}
locustag_translation = {}
locustag_COGmatch = {}
locustag_PFAMmatch = {}
locustag_TIGRmatch = {}


bioGenbankRecord = enumerate(SeqIO.parse(gbk_filename, "genbank"))
for (index, seq_record) in bioGenbankRecord:
##http://biopython.org/wiki/SeqRecord
#ncbi_bioGenbankRecord = enumerate(SeqIO.parse(ncbi_gbk_filename, "genbank"))
#for (index, seq_record) in ncbi_bioGenbankRecord:
    #print("index %i, ID = %s, length %i, with %i features"
    #      % (index, seq_record.id, len(seq_record.seq), len(seq_record.features)))
    #if len(seq_record.features) > 1:
    #    for feature in seq_record.features:
    #        pp.pprint(feature)

    ##assembly_gap
    ##source
    ##repeat_region
    ##CDS
    ##type: CDS
    ##location: [73:337](-)
    ##qualifiers:
    ##    Key: codon_start, Value: ['1']
    ##    Key: gene, Value: ['rplU_2']
    ##    Key: inference, Value: ['ab initio prediction:Prodigal:2.6', 'similar to AA sequence:UniProtKB:P26908']
    ##    Key: locus_tag, Value: ['FOLHONMB_02350']
    ##    Key: product, Value: ['50S ribosomal protein L21']
    ##    Key: transl_table, Value: ['11']
    ##    Key: translation, Value: ['MGQAIYVEKLDVEAGEKVVFDEVILVGGESTKVGAPTVAGATVEGTVEKHGKQKKVVTFKYKPKKHSHRKQGHRQPYTKVVIDAINA']
    if len(seq_record.features) > 1:
        locustag = ""
        CDS_start_pos = -1
        CDS_end_pos = -2
        for feature in seq_record.features:
          if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA' or feature.type == 'tmRNA':
            if "locus_tag" in feature.qualifiers:
                locustag = feature.qualifiers['locus_tag'][0]
                CDS_start_pos = feature.location.start.position
                CDS_end_pos = feature.location.end.position
                #print CDS_start_pos
                #print CDS_end_pos

            #151100..151093 --> 151093..151100
            #2809251..2806589 --> 2806589..2809251
            #3037667..3036778 --> 3036778..3037667
            #3385735..3385344 --> 3385344..3385735
            if "gene" in feature.qualifiers:
                #TODO: delete the last _number  : s.rfind('l')
                gene_string = feature.qualifiers['gene'][0]
                #print gene_string
                # if not find: the last char will be deleted.
                #locustag_gene[locustag]=gene_string[:gene_string.rfind('_')]
                locustag_gene[locustag]=gene_string
            if "product" in feature.qualifiers:
                locustag_product[locustag]=feature.qualifiers['product'][0]
            if "note" in feature.qualifiers:
                locustag_note[locustag]=feature.qualifiers['note'][0]
            if "translation" in feature.qualifiers:
                #print(feature.qualifiers['translation'])
                locustag_translation[locustag]=feature.qualifiers['translation'][0]
          elif (feature.type == 'COG_match') and (feature.location.start.position >= CDS_start_pos) and (feature.location.end.position <= CDS_end_pos):
            if "product" in feature.qualifiers:
                cog_string = feature.qualifiers['product'][0]
                #print(" ".join(cog_string.split(" ")[:-3]))
                #[X]COG2014: COG5527 Protein involved in initiation of plasmid replication  -3
                #[X]COG2014: COG5527 Protein involved in initiation of plasmid replication  Length=316 Score=399   -1
                #cog annotation
                locustag_COGmatch[locustag]= " ".join(cog_string.split(" ")[:-3])

               
          elif (feature.type == 'PFAM_match') and (feature.location.start.position >= CDS_start_pos) and (feature.location.end.position <= CDS_end_pos):
            if "product" in feature.qualifiers:
                if locustag in locustag_PFAMmatch:
                    locustag_PFAMmatch[locustag].append(feature.qualifiers['product'][0])
                else:
                    locustag_PFAMmatch[locustag] = [feature.qualifiers['product'][0]]
          elif (feature.type == 'TIGR_match') and (feature.location.start.position >= CDS_start_pos) and (feature.location.end.position <= CDS_end_pos):
            if "product" in feature.qualifiers:
                #print(feature.qualifiers['product'][0])
                if locustag in locustag_TIGRmatch:
                    locustag_TIGRmatch[locustag].append(feature.qualifiers['product'][0])
                else:
                    locustag_TIGRmatch[locustag] = [feature.qualifiers['product'][0]]

#pp.pprint(locustag_TIGRmatch)
#    'group_2293': [   'spoVE; GO:0030436 GO:0009252 GO:0016021 GO:0003674',
#                      'spoVE; GO:0030436 GO:0009252 GO:0016021 GO:0003674',
#                      'ftsW; GO:0051301 GO:0016021 GO:0009252 GO:0003674'],




# TODO: delete the first record!



##TODO: delete the _number at the end in Swissport_Annotation
##Swiss-Prot
## input the file gene_presence_absence.csv and output reformatted 
##"Non-unique Gene name"
lines = open(input_for_reformatting)
padj = 0.0
for line in lines:
    if not line.startswith("\"\""):
        locus_tag = line.split(",")[0].strip()
        locus_tag = locus_tag.strip("\"")    #BIT33 Gene ID
        
        #log2fc = 0.0
        #fc = 0.0
        #padj = 0.0
        log2fc = float(line.split(",")[2].strip())
        log2fc_r = round(log2fc, 2)
        fc = round(2**log2fc, 2)
        padj_str = line.split(",")[6].strip()
        if padj_str != "NA":
          padj = round(float(line.split(",")[6].strip()), 2)
        
        #print(locus_tag)
        gene_name = ""
        if locus_tag in locustag_gene:
            gene_name = locustag_gene[locus_tag]  
        product = ""
        if locus_tag in locustag_product:
            product = locustag_product[locus_tag]
        if locus_tag in locustag_note:
            note = locustag_note[locus_tag]
        translated_seq = ""
        if locus_tag in locustag_translation:
            translated_seq = locustag_translation[locus_tag]

        cogmatch = ""
        cogcode = ""
        if locus_tag in locustag_COGmatch:
            cogmatch = locustag_COGmatch[locus_tag]
            #if cogmatch:
            cogmatch_splited = cogmatch.split("]")
            cogmatch = cogmatch_splited[1]
            cogcode = cogmatch_splited[0].lstrip("[")
        pfammatch = ""
        if locus_tag in locustag_PFAMmatch:
            pfammatch = " ".join(locustag_PFAMmatch[locus_tag][0].split(" ")[:-2])
        tigrmatch = ""
        if locus_tag in locustag_TIGRmatch:
            tigrmatch = locustag_TIGRmatch[locus_tag][0].split(";")[0]

        ##(",".join(line.split(",")[3:])).strip()
        #swissprot_product, pfammatch, tigrmatch, 
        #print log2fc
        #print fc
        #print padj
        print("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%.2f,%.2f,%.2f,\"%s\""%(locus_tag, gene_name, product, note, cogcode, cogmatch, log2fc_r, fc, padj, translated_seq))
    else: 
        ##print('"GroupId","SwissProt_Annotation","NCBI_Annotation","Prokka_Annotation","SwissProt_Product","NCBI_Product","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc","HD04-03","HD4N15","Translation"')
        print('"Gene ID","Gene Name","Product","Note","COG Code","COG Annotation","log2 Fold Change","Fold Change","adj.P.Value","Translation"')
        

###      BIT33 Gene ID | Annotation | COG Annotation | log2 Fold Change | Fold Change | P-Value 
##round(2**(-0.11),2)
#print(len(locustag_COGmatch))
#print(len(locustag_PFAMmatch))
#print(len(locustag_TIGRmatch))

