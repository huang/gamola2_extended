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

## 8, an example process for an annotated scoary file
```sh
#--- 1 ----
cd /media/jhuang/Titisee/GAMOLA2
construct_Results.sh 

#---- 3.1 (complete, slow, optional) ----
cd roary_182s_95  # since the DIR Qi_panGenome contains roary results
echo ">roary_182s_95" > roary_182s_95.fa;
merge_seq.py pan_genome_reference.fa >> roary_182s_95.fa;
format_fasta_header.py pan_genome_reference.fa > pan_genome_reference_.fa
samtools faidx pan_genome_reference_.fa
bioawk -c fastx '{ print $name, length($seq) }' < pan_genome_reference_.fa > length.txt 
generate_gene_model.py length.txt > /media/jhuang/Titisee/GAMOLA2/Results/gene_models/roary_182s_95.fa.combined;
cp roary_182s_95.fa /media/jhuang/Titisee/GAMOLA2/Input_sequences
cd /media/jhuang/Titisee/GAMOLA2

#---- 3.2 (500 genes of interest, quick) ----
#select 500 genes of interest into pan_genome_reference_selected.fa
cd scoary
cp ../roary/pan_genome_reference.fa ./
format_fasta_header.py pan_genome_reference.fa > pan_genome_reference_.fa
samtools faidx pan_genome_reference_.fa
cut -f1-1 -d',' If_inf_28_07_2022_1148.results.csv > genes_of_interest.txt
#transfer the content of genes_of_interest.txt into the following gene_list
for gene in group_2034 srfAD group_2229 group_2233 group_2230 group_1857 group_1358 group_1974 group_150 group_2474 group_1360 group_2143 group_2306 group_454 group_1339 group_3239 group_622 group_1167 group_3146 group_1930 group_2103 group_2698 group_8800 group_2412 group_1756 group_2998 group_30 group_1946 melA group_1361 group_2999 group_534 ykoD_2 group_1170 oppC group_1966 group_1967 group_2399 irtB group_1538 group_1664 group_839 rl5 group_998 group_8877 group_3215 sbnB sbnA group_1536 group_1987 group_1165 group_1168 group_1535 group_1968 group_1969 glcK group_251 nanK gsiC_2 group_3120 group_455 group_621 group_850 group_1175 yxlF_2 group_877 group_1590 group_2617 group_1506 msbA_2 sufU ilvH group_1140 group_1459 group_1479 group_1512 group_1646 group_1685 atpH_2 group_1801 group_1848 sgbE prfA_1 btuF_3 group_218 group_2288 rpmE2 group_272 group_286 mnaT group_372 leuD group_500 malF malE rbfA nanT_2 yfmO group_64 group_718 dmsB cobN_2 group_807 dmsA group_9 group_1596 group_2104 dasA group_2109 group_2128 dapA_3 group_2704 group_2706 group_2707 group_2708 group_2709 group_2710 group_2712 group_2713 group_2714 group_2715 group_2716 group_2717 group_2718 group_2719 group_2721 group_2722 group_2725 lysA_1 group_2727 group_2728 group_2729 group_2734 group_8774 ulaF_1 group_8776 group_8779 folE_2 queD group_8782 group_8783 group_8784 group_8785 group_8786 group_8787 group_8788 group_8789 group_8790 group_8791 group_8793 group_8794 group_8795 group_8796 group_8797 group_8798 group_8799 group_8801 group_8802 group_8804 group_8805 group_8806 group_8807 group_8808 group_8809 group_8810 group_8811 degU group_8813 group_8814 group_8815 group_8816 group_8817 group_8818 group_8819 group_8820 tcyC group_8822 group_8823 group_8824 group_8825 group_8826 group_8827 group_8828 group_8830 group_8831 group_8832 group_8834 group_8835 group_8836 group_8837 group_8838 group_8839 group_8840 group_8841 group_8842 group_8844 glpR_4 group_8846 group_8847 dat group_8849 tdh_2 group_8851 group_8852 group_8853 group_8854 group_8855 group_8856 group_8858 group_8859 group_8860 group_8861 group_8862 group_8863 group_8864 dppC_2 group_8866 oppF_1 group_8868 group_8869 group_8870 group_8871 group_8872 group_8873 group_8874 group_8875 group_8876 group_8878 group_8879 group_8880 group_8881 group_8882 group_8883 group_8884 group_8885 btuD_9 group_8887 group_8888 group_8889 group_8890 group_8891 group_8892 group_8893 group_8895 group_8896 group_8898 group_8899 group_8901 group_8902 group_8903 group_8904 group_8905 group_8906 group_8907 group_8908 group_8909 group_8910 group_8911 group_8912 hsdR_2 hsdR_1 group_8916 btuD_10 group_8918 group_8920 group_8921 group_8922 group_8923 group_8924 group_8925 group_8926 group_8927 group_8928 group_8930 group_8931 group_8932 group_8933 group_8935 group_8936 group_8937 group_8938 group_8939 group_8940 group_8941 group_8942 group_8943 group_8952 group_8956 group_8962 group_8963 group_8964 group_8965 group_8966 group_8967 group_8968 group_8969 group_8970 yecS group_8972 group_8973 group_8974 group_8975 group_8976 group_8977 group_8978 group_8979 group_8981 group_8982 group_8983 group_8984 group_8985 leuC group_8987 group_8989 group_8990 group_8991 group_8992 group_8993 group_8994 group_8995 group_8996 group_8997 group_8998 group_8999 group_9000 group_9001 group_9002 group_9003 group_9004 group_9005 group_9006 group_9007 group_9008 group_9009 group_9010 group_9011 group_9013 group_9014 group_9015 group_9016 group_9017 group_9018 group_9019 group_9020 group_9021 group_9022 group_9023 group_9024 mbtB group_9026 group_9027 group_9028 group_9030 group_9031 yafQ group_9033 group_9411 group_1859 group_797 group_1970 group_2383 group_1164 group_1971 group_2400 group_2401 group_2232 melC_1 melD_1 group_1178 group_9519 group_2411 ortA metXA group_1213 group_2080 group_2647 group_2723 group_2724 group_1455 group_1860 menA group_407 group_1095 lolD group_3042 group_2404 group_3100 group_3126 group_851 group_852 group_882 group_2406 group_1976 iolE_1 group_2711 group_2720 group_8777 group_8833 group_8894 group_8897 group_8934 group_8944 group_8945 coaBC_1 group_8947 group_8948 group_8949 bioM_1 group_8951 group_8953 group_8954 group_8955 group_8957 group_8958 group_8959 group_8960 group_8961 group_9029 group_59 group_3265 irtA group_2314 group_3123 group_2083 group_770 group_855 group_582 group_1965 group_3065 znuC group_2231 dnaE2 group_3057 group_1978 group_1231 group_840 group_2407 group_2112 group_3119 group_1180 pdhC group_1630 metB group_3203 group_3204 group_3205 group_3206 group_3238 group_3117 group_3354 group_3406 group_8778 group_8843 group_8857 group_8919 group_9012 group_9034 group_9035 group_3102 group_571 ttuB group_1478 group_3058 group_2988 ispF group_3237 group_1319 group_1320 group_1721 group_2198 group_2989 group_2991 group_513 group_970 group_2278 group_3207 group_2133 group_2134 group_9561 group_9628 group_9698 group_224 group_293 group_844 group_1988 group_3125 group_1959 group_2444 group_3271 group_120 group_1288 group_1560 group_1713 bchI_2 yfmC group_53 dinB group_960 group_1595 group_2699 group_2700 group_2701 group_1939 group_2234 group_2357 group_3043 group_2023 group_3089 rihA group_3191 group_3751 group_2987; do
samtools faidx pan_genome_reference_.fa ${gene} >> pan_genome_reference_selected.fa;
done
echo ">roary_186_selected_genes" > roary_186_selected_genes.fa;
merge_seq.py pan_genome_reference_selected.fa >> roary_186_selected_genes.fa;
samtools faidx pan_genome_reference_selected.fa
bioawk -c fastx '{ print $name, length($seq) }' < pan_genome_reference_selected.fa > length.txt 
generate_gene_model.py length.txt > /media/jhuang/Titisee/GAMOLA2/Results/gene_models/roary_186_selected_genes.fa.combined;
cp roary_186_selected_genes.fa /media/jhuang/Titisee/GAMOLA2/Input_sequences
cd /media/jhuang/Titisee/GAMOLA2

#---- 4 ----
cd GAMOLA2
#WARNING: !!swissprot as default database, manually choose nr as Blast_db --> too slow, better choose swissprot, it is quick!!
./Gamola.pl    #No Glimmer model and Critica database due to self-extracted ORF; choosing nr.pal or swissprot.pal as Blast_db; COG2014 as COG_db; Pfam-A.hmm as Pfam_db; No Rfam_db; TIGRFAMS_15.0_HMM.LIB as TIGRfam_db
/media/jhuang/Titisee/GAMOLA2/Results/COG_results
#grep "Class: ," roary_182s_95.fa_COG_*
for sample in 10601 10710 11187 11377 12474 12504  1520 1551 2092 2160 2467   289 3455 4694 5862 5863 6166 6464 7618 8007 8406 8784 9157 9373 953; do
  sed -i -e 's/Class: ,/Class: None, None/g' roary_182s_95.fa_COG_${sample}
done
./Gamola.pl

#---- 5 ---- update locustag
cp /media/jhuang/Titisee/GAMOLA2/Consolidated_results/roary_182s_95.fa/roary_182s_95.fa.gb ./
awk '{print $0 "ORIGIN"> "file" NR}' RS='ORIGIN' roary_182s_95.fa.gb
sed -i 's/</ /g' file2
sed -i 's/^ //g' file2
cat file1 file2 > roary_182s_95_.fa.gb
iconv -t UTF-8 -f Windows-1252 roary_182s_95_.fa.gb
#89600270 Jan 26 13:44 roary_182s_95__.fa.gb
#-- if ran [3.1 complete, optional]--
python ~/Scripts/update_locustag.py roary_182s_95_c.fa.gb /home/jhuang/DATA/Data_Anna_C.acnes/182samples_roaries/roary_182s_95/pan_genome_reference_.fa.fai 
#-- if ran [3.2 selected]--
python ~/Scripts/update_locustag.py roary_186_selected_genes_.fa.gb /home/jhuang/DATA/Data_Anna_C.acnes/scoary/pan_genome_reference_selected.fa.fai > roary_186_selected_genes__.fa.gb 
> roary_186_selected_genes__.fa.gb    #clean old *__.fa.gb
#DEBUG: repeated record "exo" in pan_genome_reference_.fa.fai, should be renamed!

#---- 6.6 ----
process_gamola_gb_scoary.py If_inf_02_09_2021_1151.results.csv ../../gamola2/roary_182s_95__.fa.gb > annotated_scoary.txt
Some fields may be wrong.
  BiopythonParserWarning)
/usr/local/lib/python2.7/dist-packages/Bio/GenBank/__init__.py:1047: BiopythonParserWarning: Ignoring invalid location: '1186620..1185596'
sed -i -e 's/1186620\.\.1185596/1185596\.\.1186620/g' roary__.fa.gb
#DEBUG "group_12676" -> "group_12667"

#---- ADD DNA-Sequences add the end of table as the last column ----
cut -d',' -f1-1 annotated_scoary.txt > get_seq_ORF.sh
#extend the file get_seq_ORF.sh to extract all DNA-sequences
#mv the bash file under DIR of pan_genome_reference_.fa.fai 
#samtools faidx pan_genome_reference_.fa group_2201 > seq_ORF.fasta
#samtools faidx pan_genome_reference_.fa group_2205 >> seq_ORF.fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < seq_ORF.fasta > seq_ORF_.fasta
grep -v ">" seq_ORF_.fasta > seq_ORF__.txt
paste -d',' annotated_scoary.txt seq_ORF__.txt > annotated_scoary_.txt
```
