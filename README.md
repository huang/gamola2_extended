# extending gamola2

#HD04/seq-seq-pan_prokka/ssp_cons.gbk


## 1, construct env under env bengal2
```sh
#under /media/jhuang/Titisee/GAMOLA2
construct_Results.sh 
#clean_Results.sh
```


## 2, assembly and prokka two assembly.fa
```sh
# under bengal3_ac3
ln -s ~/REFs/plasmid_databases databases
grep "Staphylococcus" plasmids.csv > plasmids_Staphylococcus.csv  # add the header-line
hyasp create ncbi_Staphylococcus_genes.fasta -p databases/plasmids_Enterococcus.csv -b databases/ncbi_blacklist.txt -d -l 500 -m 100 -t GenBank
hyasp_pipeline HJB-2-1_hyasp_out ncbi_Staphylococcus_genes.fasta -1 ~/DATA/Data_Anna13_vanA-Element/raw_data/HJB-2-1/642_01_S23_R1_001.fastq -2 ~/DATA/Data_Anna13_vanA-Element/raw_data/HJB-2-1/642_01_S23_R2_001.fastq    # NOTE no fastq.gz input

ln -s ~/REFs/plasmid_databases databases
hyasp create ncbi_Staphylococcus_genes.fasta -p databases/plasmids_Enterococcus.csv -b databases/ncbi_blacklist.txt -d -l 500 -m 100 -t GenBank
for sample in HD04_10_S79 HD04_1_S70 HD04_2_S71 HD04_3_S72 HD04_4_S73 HD04_5_S74 HD04_6_S75 HD04_7_S76 HD04_8_S77 HD04_9_S78 HD12_10_S89 HD12_1_S80 HD12_2_S81 HD12_3_S82 HD12_4_S83 HD12_5_S84 HD12_6_S85 HD12_7_S86 HD12_8_S87 HD12_9_S88 HD15_10_S10 HD15_1_S1 HD15_2_S2 HD15_3_S3 HD15_4_S4 HD15_5_S5 HD15_6_S6 HD15_7_S7 HD15_8_S8 HD15_9_S9 HD17_1_S90 HD17_2_S91 HD17_3_S92 HD17_4_S93 HD17_5_S94 HD17_6_S95 HD17_7_S96 HD21N14_S67 HD25_11_S11 HD25_12_S12 HD26_1_S69 HD40_5_S16 HD46_8_S68 HD4N1_S13 HD4N2_S14 HD4N6_S15 HD59_10_S66 HD59_1_S57 HD59_2_S58 HD59_3_S59 HD59_4_S60 HD59_5_S61 HD59_6_S62 HD59_7_S63 HD59_8_S64 HD59_9_S65; do
hyasp_pipeline 180119_N_hyasp/${sample} ncbi_Staphylococcus_genes.fasta -1 /media/jhuang/Elements/Data_Anna12_HAPDICS_HyAsP/180119_N/${sample}_R1_001.fastq -2 /media/jhuang/Elements/Data_Anna12_HAPDICS_HyAsP/180119_N/${sample}_R2_001.fastq
done

#--racon --keep_subplasmids 
for sample in HD04-03 HD4N15 HD04-08 HD4N01 HD21-7 HD21N14 HD26-5 HD26N1 HD27-3 HD27N1 HD29-5 HD29N4 HD33-2 HD33N22 HD59-01 HD59N26; do
prokka --force --outdir prokka --cpus 10 --kingdom Bacteria --genus Staphylococcus --species epidermidis --addgenes --addmrna --prefix ${sample} --locustag ${sample} ${sample}_hyasp_out/assembly/assembly.fasta
done
roary -p 14 -f ./roary -i 95 -cd 99 -s -e -n -v prokka/HD04-03.gff prokka/HD4N15.gff
mv roary roary_HD04_comp1
roary -p 14 -f ./roary -i 95 -cd 99 -s -e -n -v prokka/HD04-08.gff prokka/HD4N01.gff
mv roary roary_HD04_comp2
roary -p 14 -f ./roary -i 95 -cd 99 -s -e -n -v prokka/HD21-7.gff prokka/HD21N14.gff 
mv roary roary_HD21_comp
roary -p 14 -f ./roary -i 95 -cd 99 -s -e -n -v prokka/HD26-5.gff prokka/HD26N1.gff
mv roary roary_HD26_comp
roary -p 14 -f ./roary -i 95 -cd 99 -s -e -n -v prokka/HD27-3.gff prokka/HD27N1.gff
mv roary roary_HD27_comp
roary -p 14 -f ./roary -i 95 -cd 99 -s -e -n -v prokka/HD29-5.gff prokka/HD29N4.gff 
mv roary roary_HD29_comp
roary -p 14 -f ./roary -i 95 -cd 99 -s -e -n -v prokka/HD33-2.gff prokka/HD33N22.gff
mv roary roary_HD33_comp
roary -p 14 -f ./roary -i 95 -cd 99 -s -e -n -v prokka/HD59-01.gff prokka/HD59N26.gff
mv roary roary_HD59_comp
# move them to HD04_comp1/prokka and roary DIR.

#prokka --usegenus --kingdom Bacteria --genus Staphylococcus --species epidermidis --strain ssp_cons --addgenes --addmrna --outdir seq-seq-pan_prokka --prefix ssp_cons --locustag ssp_cons --force seq-seq-pan/seq-seq-pan_consensus.fasta.blockseparated.fasta
#-hmm /media/jhuang/Titisee/GAMOLA2/TIGRfam_db/TIGRFAMs_15.0_HMM.LIB
```


## 3.1, preparing gene_models from roary
```sh
for sample in HD04_comp1 HD04_comp2 HD21_comp HD26_comp HD27_comp HD29_comp HD33_comp HD59_comp; do
    cd ${sample}/roary
    echo ">${sample}" > ${sample}.fa;
    merge_seq.py pan_genome_reference.fa >> ${sample}.fa;
    format_fasta_header.py pan_genome_reference.fa > pan_genome_reference_.fa
    samtools faidx pan_genome_reference_.fa
    bioawk -c fastx '{ print $name, length($seq) }' < pan_genome_reference_.fa > length.txt 
    generate_gene_model.py length.txt > ../../GAMOLA2/Results/gene_models/${sample}.fa.combined;
    cp ${sample}.fa ../../GAMOLA2/Input_sequences
    cd ../..
done

#under ALL78
for sample in roary; do
    cd ${sample}  # since the DIR Qi_panGenome contains roary results
    echo ">${sample}" > ${sample}.fa;
    merge_seq.py pan_genome_reference.fa >> ${sample}.fa;
    format_fasta_header.py pan_genome_reference.fa > pan_genome_reference_.fa
    samtools faidx pan_genome_reference_.fa
    bioawk -c fastx '{ print $name, length($seq) }' < pan_genome_reference_.fa > length.txt 
    generate_gene_model.py length.txt > /media/jhuang/Titisee/GAMOLA2/Results/gene_models/${sample}.fa.combined;
    cp ${sample}.fa /media/jhuang/Titisee/GAMOLA2/Input_sequences
    cd ..
done
```

## 3.2(deprecated), preparing gene_models from seq-seq-pan_prokka/ssp_cons.ffn
```sh
for sample in HD04 HD05 HD104 HD12 HD15 HD17 HD21 HD25 HD26 HD27 HD29 HD31 HD33 HD39 HD40 HD43 HD46 HD47 HD59 HD66 HD69 HD75 HD99; do
for sample in HD33; do
    cd ${sample}/seq-seq-pan_prokka
    echo ">${sample}" > ${sample}.fa;
    merge_seq.py ssp_cons.ffn >> ${sample}.fa;
    format_fasta_header_using_part0.py ssp_cons.ffn > ssp_cons_.fa
    samtools faidx ssp_cons_.fa
    bioawk -c fastx '{ print $name, length($seq) }' < ssp_cons_.fa > length.txt 
    generate_gene_model.py length.txt > ../../GAMOLA2/Results/gene_models/${sample}.fa.combined;
    cp ${sample}.fa ../../GAMOLA2/Input_sequences
    cd ../..
done
```

## 3.3(doesn't work!), preparing gene_models from INPUT_FILE[seq-seq-pan_prokka] for variants_clonal.xls (from gff-file to fa.combined)
```sh
# #Assigned 4964 locus_tags to CDS and RNA features (220 tRNAs + 7 rRNAs + 4737 CDS = 227 + 4737 = 4964), and + 3 CRISPRs
# for sample in HD04 HD05 HD104 HD12 HD15 HD17 HD21 HD25 HD26 HD27 HD29 HD31 HD33 HD39 HD40 HD43 HD46 HD47 HD59 HD66 HD69 HD75 HD99; do
# for sample in HD33; do
# cd ${sample}/seq-seq-pan_prokka
# grep -P "seq-seq-pan\tprokka\tgene" ssp_cons.gff > ssp_cons_gene.gff
# generate_gene_model_from_gff.py ssp_cons_gene.gff > ../../GAMOLA2/Results/gene_models/${sample}.fa.combined
# cp ssp_cons.fsa ../../GAMOLA2/Input_sequences/${sample}.fa
# cd ../../
# done

# cd GAMOLA2/Input_sequences
# for sample in HD04 HD05 HD104 HD12 HD15 HD17 HD21 HD25 HD26 HD27 HD29 HD31 HD33 HD39 HD40 HD43 HD46 HD47 HD59 HD66 HD69 HD75 HD99; do
# for sample in HD33; do
# # under gamola2/Results/gene_models
# sed -i -e 's/ssp_cons_0000//g' ${sample}.fa.combined
# sed -i -e 's/ssp_cons_000//g' ${sample}.fa.combined
# sed -i -e 's/ssp_cons_00//g' ${sample}.fa.combined
# sed -i -e 's/ssp_cons_0//g' ${sample}.fa.combined
# done
```

## 3.4, preparing gene_models for RNASeq (input RP62A.gb, RP62A.fa, RP62A.gff → RP62A_cds.fa and RP62A_cds.model)
```sh
#gamola2 with RNASeq-output: e.g. under /media/jhuang/Elements/Data_Tam_RNASeq
# under results/reference_genome 
#paaI 
#https://www.biostars.org/p/222028/
#-- fasta -- 
grep "ID=gene" RP62A.gff > RP62A_.gff 
bedtools getfasta -fi RP62A.fa -bed RP62A_.gff -s -fo RP62A_.fna -name 
cut -f1-2 -d$'\t' RP62A_.gff > f1-2 
cut -f4- -d$'\t' RP62A_.gff > f4- 
cut -f9-9 -d$'\t' RP62A_.gff > f9 
cut -f2-2 -d';' f9 > f9_2 
paste -d$"\t" f1-2 f9_2 f4- > RP62A__.gff 
grep -v "pseudo=true" RP62A__.gff > RP62A___.gff 
#grep -v "gene_biotype=RNase_P_RNA" 
grep -v "RNA" RP62A___.gff > RP62A____.gff   #2528 vs 2536
#oder depending on the input gff format
#grep "gene_biotype=protein_coding" RP62A.gff > RP62A____.gff   #2526
bedtools getfasta -fi RP62A.fa -bed RP62A____.gff -s -fo RP62A.fna -name
sed -i -e 's/Name=//g'  RP62A.fna

echo ">RP62A_cds" > RP62A_cds.fa; 
merge_seq.py RP62A.fna >> RP62A_cds.fa; 

#-- fna.fai – 
samtools faidx RP62A.fna
#NOTE that samtools faidx GCF_002811175.fna;  GCF_002811175.fna.fai has some bugs --> not correct!!!

#-- gene_model --
bioawk -c fastx '{ print $name, length($seq) }' < RP62A.fna > length.txt 
#check the length are always 3 times of nt, [BIT33_RS12745   1096] is not correct! 
# search delete the first base in BIT33_RS12745 in GCF_002811175.fna: rerunning this block. 
# WKLILI*TN*KT* LIVAKHYGGIFDYDLKK,  DELETE 39 nt bases at the beginning in GCF_002811175.fna. 1056
#comparing the protein with #https://www.ncbi.nlm.nih.gov/protein/WP_071543072.1
report=genbank&log$=prottop&blast_rank=1&RID=H1772H02014 

generate_gene_model.py length.txt > GAMOLA2/Results/gene_models/RP62A_cds.fa.combined;
#generate_gene_model.py length.txt > model.temp4 
#cp model.temp4 ../../GAMOLA2/Results/gene_models/BIT33.fa.combined; 
cp RP62A_cds.fa GAMOLA2/Input_sequences

# MODIFYING the start_codons in source code 
#grep "start" Error.log | sort -u > no_classic_startcodon.txt 
#ata, atc, att, ctg
```


## 4, running Gamola.pl
```sh
#           for HD04, HD05, ...: running GAMOLA2, select only one type annotation, for example blastp (Swissprot or NCBI) --> HD04_comp1_SwissProt__.fa.gb and HD04_comp1_NCBI__.fa.gb
#https://www.biostars.org/p/336016/
cd GAMOLA2
# IMPORTANT: set the nr_Staphy as the Blast_db
./Gamola.pl
# add more tolerated start and top codons into lib/ProgrammeModules/sequence.pm

#BUG: Error parsing COG2014 results for file HD04_comp1.fa, see GAMOLA2/Error.log
Could not parse the COG2014 code from entry YP_002559933    COG2014: COG5604; Class: ,                                                               ;                   ; Annotation:  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Length=2045  Score=1196  Expect=0.0.
http://prodata.swmed.edu/citrusgreening/psyllid/cogdata/COG5604_genelist.html
#SOLUTION under GAMOLA2/Results: sed -i -e 's/Class: ,/Class: None, None/g' Results/COG_results/roary_80.fa_COG_*
#HD04_comp1.fa_COG_1022:>YP_002559933    COG2014: COG5604; Class: None, None                                                   ; Genome: Macrococcus_caseolyticus_JCSC5402                ; Annotation:  None~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Length=2045  Score=3  Expect=0.015

#under GAMOLA2/Results/COG_results
grep "Class: ," HD04.fa_COG*
for sample in 1489 1833 2050 2236 2309 321 3282 3306 3491 3693 4903 4904 4906 556 574 611 637 638 644 648 712 935; do
sed -i -e 's/Class: ,/Class: None, None/g' HD04.fa_COG_${sample}
done
grep "Class: ," HD05.fa_COG*
for sample in 1206 137 1371 1425 1433 1465 1603 2249 229 2320 2569 557 833; do
sed -i -e 's/Class: ,/Class: None, None/g' HD05.fa_COG_${sample}
done
grep "Class: ," HD12.fa_COG*
for sample in 1158 1170 126 1830 2002 2045 2289 2947 328 466 501 520 577 699; do
sed -i -e 's/Class: ,/Class: None, None/g' HD12.fa_COG_${sample}
done

grep "Class: ," HD15.fa_COG*
for sample in 1203 13 167 1852 2284 2538 2710 3032 3037 608 705 740 798 802 824; do
sed -i -e 's/Class: ,/Class: None, None/g' HD15.fa_COG_${sample}
done
grep "Class: ," HD17.fa_COG*
for sample in 1763 1925 3235 4535 5702 963; do
sed -i -e 's/Class: ,/Class: None, None/g' HD17.fa_COG_${sample}
done

grep "Class: ," HD21.fa_COG*
for sample in 1542 159 1757 1791 1794 1797 1910 1919 2426 2743 2920 3226 3846 3912 58; do
sed -i -e 's/Class: ,/Class: None, None/g' HD21.fa_COG_${sample}
done

grep "Class: ," HD25.fa_COG*
for sample in 1498 1511 1569 1572 1853 1881 2015 2095 2810 3344 3547 3780 422 519 827 986; do
sed -i -e 's/Class: ,/Class: None, None/g' HD25.fa_COG_${sample}
done
grep "Class: ," HD26.fa_COG*
for sample in 1211 1257 1340 1342 1350 1528 1908 2400 3766 3775 3804 4023 4030 4376 4483 695 823; do
sed -i -e 's/Class: ,/Class: None, None/g' HD26.fa_COG_${sample}
done
grep "Class: ," HD27.fa_COG*
for sample in 1011 1073 1529 1532 1569 1639 1645 2571 2867 2869 3384 3668 3729 3731 4241 4462 4501 545 625; do
sed -i -e 's/Class: ,/Class: None, None/g' HD27.fa_COG_${sample}
done
grep "Class: ," HD29.fa_COG*
for sample in 1007 148 1663 1915 196 2193 24 2569 2731 617 843 849 852 908 944; do
sed -i -e 's/Class: ,/Class: None, None/g' HD29.fa_COG_${sample}
done
grep "Class: ," HD31.fa_COG*
for sample in 1648 1723 1776 2001 205 2129 2148 2357 2371 2622 460 638 70 807; do
sed -i -e 's/Class: ,/Class: None, None/g' HD31.fa_COG_${sample}
done
grep "Class: ," HD39.fa_COG*
for sample in 1010 1487 1732 2326 257 2586 2991 412 509 512 619 768; do
sed -i -e 's/Class: ,/Class: None, None/g' HD39.fa_COG_${sample}
done
grep "Class: ," HD40.fa_COG*
for sample in 1296 1631 1738 1847 1883 2055 2477 2739 276 2933 330 639 697 700; do
sed -i -e 's/Class: ,/Class: None, None/g' HD40.fa_COG_${sample}
done
grep "Class: ," HD43.fa_COG*
for sample in 100 1405 1689 1732 1800 1804 2354 2451 2453 2454 2469 2561 2562 2563 3012 3099 3100 3103 3104 3105 3164 3557 3580 3622 3714 3742 3753 3925 449 589 931; do
sed -i -e 's/Class: ,/Class: None, None/g' HD43.fa_COG_${sample}
done
grep "Class: ," HD46.fa_COG*
for sample in 141 1572 2119 234 2437 2578 2752 2900 292 450 458 512 545 561 738; do
sed -i -e 's/Class: ,/Class: None, None/g' HD46.fa_COG_${sample}
done

grep "Class: ," HD47.fa_COG*
for sample in 1365 1899 2008 208 2382 2728 2767 3146 3269 736 79 834 837 893 927; do
sed -i -e 's/Class: ,/Class: None, None/g' HD47.fa_COG_${sample}
done
grep "Class: ," HD59.fa_COG*
for sample in 1001 1004 1005 147 1681 1786 2022 2042 206 2944 2990 584 825 901 938; do
sed -i -e 's/Class: ,/Class: None, None/g' HD59.fa_COG_${sample}
done
grep "Class: ," HD66.fa_COG*
for sample in 1043 1199 1202 1261 1713 1826 2282 2434 2483 2718 2890 2910 2911 2953 2954 2959 2971 2976 3028 3029 3051 3058 3063 3070 3077 3078 3087 3164 3240 3242 3501 3713 411 511 955 970; do
sed -i -e 's/Class: ,/Class: None, None/g' HD66.fa_COG_${sample}
done
grep "Class: ," HD69.fa_COG*
for sample in 1013 1059 1060 1062 1164 1214 1245 1737 1942 1981 2236 2653 2654 2793 2846 2847 439 499; do
sed -i -e 's/Class: ,/Class: None, None/g' HD69.fa_COG_${sample}
done
grep "Class: ," HD75.fa_COG*
for sample in 1067 1745 1748 1812 1846 1862 1904 2269 2305 2420 2421 3091 3092 3093 313 3185 3267 3268 3269 3316 3510 3512 3585 76; do
sed -i -e 's/Class: ,/Class: None, None/g' HD75.fa_COG_${sample}
done
grep "Class: ," HD99.fa_COG*
for sample in 1192 141 1563 1612 1836 1854 2290 235 2396 2465 358 643 717 738 790 823 824 90; do
sed -i -e 's/Class: ,/Class: None, None/g' HD99.fa_COG_${sample}
done
grep "Class: ," HD104.fa_COG*
for sample in 1421 1632 3124 4221 4222 4246 457 4833 4909 4929 5541 5717 5722 5813 5846 6294 7372 7384 7597 7701 7903 7933 7934 8349 8545 9516 9619; do
sed -i -e 's/Class: ,/Class: None, None/g' HD104.fa_COG_${sample}
done
grep "Class: ," HD33.fa_COG*  # 102    22K Sep 18 11:24 HD33.fa
for sample in 1099 1103 1170 1604 1817 1862 1928 220 2366 2595 3280 3493 3694 3719 4013 4154 4204 4211 428 4291 4324 4499 4523 4542 4558 4561 4623 549; do
sed -i -e 's/Class: ,/Class: None, None/g' HD33.fa_COG_${sample}
done

./Gamola.pl
```


## 5, update_locustag.py under e.g. ~/DATA/Data_Anna12_HAPDICS_final/
```sh
#for HD**: 
for sample in HD04 HD05 HD104 HD12 HD15   HD17 HD21 HD25 HD26 HD27   HD29 HD31 HD33 HD39 HD40   HD43 HD46 HD47 HD59 HD66   HD69 HD75 HD99; do
#for HD**_comp and roary: 
for sample in roary; do
#for sample in RP62A_cds; do   #needs ‘mkdir RP62A_cds’ under featureCounts having DIRs GAMOLA2 and degenes
# under the DIR ./gamola2, update_locustag using "../roary/pan_genome_reference_.fa.fai" or "../seq-seq-pan_prokka/ssp_cons_.fa.fai" generated in 3.1.1
mkdir ${sample}/gamola2;
cp /media/jhuang/Titisee/GAMOLA2/Consolidated_results/${sample}.fa/${sample}.fa.gb ./${sample}/gamola2;
cd ${sample}/gamola2;
awk '{print $0 "ORIGIN"> "file" NR}' RS='ORIGIN' ${sample}.fa.gb
sed -i 's/</ /g' file2
sed -i 's/^ //g' file2
cat file1 file2 > ${sample}_.fa.gb
#update_locustag.py ${sample}_.fa.gb ../seq-seq-pan_prokka/ssp_cons_.fa.fai > ${sample}__.fa.gb
#update_locustag.py ${sample}_.fa.gb ../../../../RP62A.fna.fai > ${sample}__.fa.gb
update_locustag.py ${sample}_.fa.gb ../../${sample}/pan_genome_reference_.fa.fai > ${sample}__.fa.gb
#NOTE that __.fa.gb suffices for further usage!
# following step retains only the feature.type = [“CDS”]
#NOT_REQUIRED: python ~/Scripts/parseGenbank_extractGenes.py -g ${sample}__.fa.gb -o ${sample}___.fa.gb
# following step corrects the genbank-format. Note that some COG_match records contain false start..end, e.g. COG_match    287381..286861 → needs to be correctly manually
#NOT_REQUIRED: python ~/Scripts/parseGenbank_reformat.py -g ${sample}__.fa.gb -o ${sample}____.fa.gb
cd ../..
done
## manually clean the empty lines after 'ORIGIN' and the last line containing 'ORIGIN' in ${sample}/gamola2/${sample}__.fa.gb
```


## 6.1, generating the annotation file with HD**_comp
```sh
# Input-format: roary/gene_presence_absence.csv
# Output-format: "Group ID","NCBI Annotation","Root","COG Code","COG Annotation","PFAM Annotation","Prokka Annotation","No. isolates","Isolate1","Isolate2","Translation"'
# CORRECT RESULTS, checking the sample order in ../roary/gene_presence_absence.csv, e.g. "HD26-5","HD26N1"
process_gamola_gb_HAPDICS.py gene_presence_absence.csv HD04_comp1__.fa.gb
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD04_comp1__.fa.gb ../prokka/HD04-03/HD04-03.gff ../prokka/HD4N15/HD4N15.gff > annotated_gene_presence_absence.csv
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD04_comp2__.fa.gb ../prokka/HD04-08/HD04-08.gff ../prokka/HD4N01/HD4N01.gff > annotated_gene_presence_absence.csv
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD21_comp__.fa.gb ../prokka/HD21-7/HD21-7.gff ../prokka/HD21N14/HD21N14.gff > annotated_gene_presence_absence.csv
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD26_comp__.fa.gb ../prokka/HD26-5/HD26-5.gff ../prokka/HD26N1/HD26N1.gff > annotated_gene_presence_absence.csv
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD27_comp__.fa.gb ../prokka/HD27-3/HD27-3.gff ../prokka/HD27N1/HD27N1.gff > annotated_gene_presence_absence.csv
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD29_comp__.fa.gb ../prokka/HD29-5/HD29-5.gff ../prokka/HD29N4/HD29N4.gff > annotated_gene_presence_absence.csv
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD33_comp__.fa.gb ../prokka/HD33-2/HD33-2.gff ../prokka/HD33N22/HD33N22.gff > annotated_gene_presence_absence.csv
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD59_comp__.fa.gb ../prokka/HD59-01/HD59-01.gff ../prokka/HD59N26/HD59N26.gff > annotated_gene_presence_absence.csv
process_gamola_gb_HAPDICS_plus_contigid.py ../roary/gene_presence_absence.csv HD59_comp__.fa.gb ../prokka/HD59-01/HD59-01.gff ../prokka/HD59N26/HD59N26.gff > annotated_gene_presence_absence.csv
# with libreoffice 'save as csv' can delete all sign '"'
```

## 6.2, generating the annotation file with HD**
```sh
# To security of the variants are based on the current LOCUS_TAG (z.g. ssp_cons_01348). Update all snippy.core.tab__ by using "rerun_clonal.sh"
# Input-format: ../variants_clonal/snippy.core.tab__
sed -i -e 's/>ssp_cons_/ssp_cons_/g'  HD99__.fa.gb
process_gamola_gb_HAPDICS_plus_variants.py ../variants_clonal/snippy.core.tab__ HD99__.fa.gb > annotated_gene_variants.csv
```

## 6.3, generating the annotation file with DESeq output
```sh
# Input-format: Data_Tam_RNASeq/format("","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj") generated by R-scripts
# Output-format: "BIT33 Gene ID","Swiss-Prot Annotation","COG Code","COG Annotation","PFAM Annotation","log2 Fold Change","Fold Change","adj.P.Value","Translation"
process_gamola_gb.py adeIJ_ko_vs_wt_background.txt BIT33__.fa.gb > annotated_degenes.csv
process_gamola_gb_raw.py d_raw.csv results/gamola2/acb_all__.fa.gb > annotated_d_raw.csv
# under bengal3_ac3
mkdir degenes_annotated;
for comp in HD04.TSB_vs_HD04.hS HD04N.hS_vs_HD04.hS HD04N.TSB_vs_HD04.TSB HD04N.TSB_vs_HD04N.hS HD04N6.hS_vs_HD04N2.hS HD04N6.TSB_vs_HD04N2.TSB HD04_26.TSB_vs_HD04_26.hS HD04_8.hS_vs_HD04_1.hS HD21N.TSB_vs_HD21.TSB HD21N.TSB_vs_HD21N.hS HD21N5.hS_vs_HD21N4.hS HD26N.hS_vs_HD26.hS HD26N.TSB_vs_HD26.TSB HD26N.TSB_vs_HD26N.hS HD04_8.TSB_vs_HD04_1.TSB HD21.TSB_vs_HD21.hS HD21N.hS_vs_HD21.hS HD21N5.TSB_vs_HD21N4.TSB HD21_2.TSB_vs_HD21_1.TSB HD26.TSB_vs_HD26.hS HD26N3.hS_vs_HD26N1.hS HD26N3.TSB_vs_HD26N1.TSB HD26_7.hS_vs_HD26_1.hS HD26_7.TSB_vs_HD26_1.TSB nose.hS_vs_infection.hS nose.TSB_vs_infection.TSB TSB_vs_hS; do \
mkdir degenes_annotated/${comp}_output
process_gamola_gb.py degenes/${comp}_output/upregulated_filtered RP62A_cds/gamola2/RP62A_cds__.fa.gb > degenes_annotated/${comp}_output/upregulated_filtered
process_gamola_gb.py degenes/${comp}_output/downregulated_filtered RP62A_cds/gamola2/RP62A_cds__.fa.gb > degenes_annotated/${comp}_output/downregulated_filtered
process_gamola_gb.py degenes/${comp}_output/background RP62A_cds/gamola2/RP62A_cds__.fa.gb > degenes_annotated/${comp}_output/background
cd degenes_annotated/${comp}_output;
~/Tools/csv2xls-0.4/csv_to_xls.py upregulated_filtered downregulated_filtered background -d$',' -o ../${comp}_annotated_degenes.xls;
cd ../..;
done
#delete the empty three columns.
```

## 6.4, generating the annotation file with roary output
```sh
# Input-format: roary/gene_presence_absence.csv
# Output-format: "Group ID","NCBI Annotation","Root","COG Code","COG Annotation","PFAM Annotation","Prokka Annotation","No. isolates","Isolate1","Isolate2","Translation"'
# CORRECT RESULTS
process_gamola_gb_roary.py ../roary/gene_presence_absence.csv HD04_comp1__.fa.gb > annotated_gene_presence_absence.csv
# Output-format: "Group ID","NCBI Annotation","Root","COG Code","COG Annotation","Prokka Annotation","Translation"'
process_gamola_gb_roary2.py gene_presence_absence.csv gamola2/roary__.fa.gb > annotated_gene_presence_absence.csv
process_gamola_gb_scoary.py If_inf_05_05_2020_1757.results.csv roary/gamola2/roary__.fa.gb > annotated_If_inf_05_05_2020_1757.results.csv
# with libreoffice 'save as csv' can delete all sign '"'
```

## 6.5, generating the annotation file with mergedSNP output
```sh
# Input-format: variants/merged_SNPs_.tab
# Output-format: CHR     POS     HD5-1   HD5-10  HD5-2   HD5-3   HD5-4   HD5-5   HD5-6   HD5-7   HD5-8   HD5-9   LOCUS_TAG       GENE    PRODUCT Effect  Impact  Functional_Class        Codon_change    Amino_Acid_change       Amino_Acid_Length       COG Code        COG Annotation  PFAM Annotation TIGR Annotation Translation
# Command:
process_gamola_gb_HAPDICS_plus_variants.py merged_SNPs_.tab ../gamola2/HD05__.fa.gb > merged_SNPs__.tab

#BUG:
'589777..589257'
  BiopythonParserWarning)
/home/jhuang/anaconda3/lib/python2.7/site-packages/Bio/GenBank/__init__.py:1047: BiopythonParserWarning: Ignoring invalid location: '753716..751369'
  BiopythonParserWarning)
/home/jhuang/anaconda3/lib/python2.7/site-packages/Bio/GenBank/__init__.py:1047: BiopythonParserWarning: Ignoring invalid location: '909573..909560'
  BiopythonParserWarning)
/home/jhuang/anaconda3/lib/python2.7/site-packages/Bio/GenBank/__init__.py:1047: BiopythonParserWarning: Ignoring invalid location: '1640188..1639431'
```

## 6.6, generating the annotation of scoary output
```sh
process_gamola_gb_scoary.py If_inf_16_08_2021_1634.results.csv roary__.fa.gb > annotated_scoary.txt
Some fields may be wrong.
  BiopythonParserWarning)
/usr/local/lib/python2.7/dist-packages/Bio/GenBank/__init__.py:1047: BiopythonParserWarning: Ignoring invalid location: '1186620..1185596'
sed -i -e 's/1186620\.\.1185596/1185596\.\.1186620/g' roary__.fa.gb
```

## 7.1, csv to xls
```sh
for sample in HD04 HD05 HD12 HD15 HD17 HD21 HD25 HD26 HD27 HD29 HD31 HD33 HD39 HD40 HD43 HD46 HD47 HD59 HD66 HD69 HD75 HD99 HD104; do
  cp ${sample}/gamola2/annotated_gene_variants.csv annotated_variants_clonal/${sample}.csv
done
#under bengal3_ac3
~/Tools/csv2xls-0.4/csv_to_xls.py HD04.csv HD05.csv HD12.csv HD15.csv HD17.csv HD21.csv HD25.csv HD26.csv HD27.csv HD29.csv HD31.csv HD33.csv HD39.csv HD40.csv HD43.csv HD46.csv HD47.csv HD59.csv HD66.csv HD69.csv HD75.csv HD99.csv HD104.csv -d$'\t' -o annotated_variants_clonal.xls
#DEBUG: 1-(5-phosphoribosyl)-5-
#group_596,1-(5-phosphoribosyl)-5-[(5-phosphoribosylamino)methylideneamino] imidazole-4-carboxamide isomerase ,Staphylococcus,E,COG2014: COG0106 Phosphoribosylformimino-5- aminoimidazole carboxamide ribonucleotide (ProFAR) isomerase ,"His_biosynth, Histidine biosynthesis protein",TIGR00007: 1-(5-phosphoribosyl)-5-[(5-phosphoribosylamino)methylideneamino]imidazole-4-carboxamide isomerase,2,HD04-03_01033,HD4N15_01036,MIDLWPAIDLINSTSVRLTEGKYDTKEKMEKSVEDSIRFYSQFKCVKRIHIVDLIGAKAKEVKEFDYIRSLRKVTTKPIEVGGGIRSKQTIENYIHSGIDYCIVGTKGIQDIEWLTHMTHQFPNKLYLSVDAFGEKIKINGWKEDAKLNLFDYVAKIEHLPLGGVIYTDISKDGKLSGPNFDLTGRLALYTSLPVIASGGIRHQEDLFRLESLNVHAAIVGKAAHLDEFWEGLS*

for sample in HD04_comp1 HD04_comp2 HD21_comp HD26_comp HD27_comp HD29_comp HD33_comp HD59_comp; do
  cp ${sample}/gamola2/annotated_gene_presence_absence_.csv annotated_gene_clusters/${sample}.csv
done
#under bengal3_ac3
~/Tools/csv2xls-0.4/csv_to_xls.py HD04_comp1.csv HD04_comp2.csv HD21_comp.csv HD26_comp.csv HD27_comp.csv HD29_comp.csv HD33_comp.csv HD59_comp.csv -d$'\t' -o annotated_gene_clusters.xls
```

## 7.2
```sh
##using libreoffice transfer the separator to '\t' and save as "annotated_gene_presence_absence_.csv"
```
