CC=gcc -Wall -O3 
LIBS= 

lookat=lookat.o Sequence.o
motiffind=motiffind.o Sequence.o
translate=translate.o Sequence.o GeneticCode.o
extractcoding=extractcoding.o Sequence.o GeneticCode.o 
intergenic=intergenic.c Sequence.o SeqMap.o
empirical-matrix=empirical-matrix.o Sequence.o GeneticCode.o
dicodontable=dicodontable.o Sequence.o SeqMap.o
addlongorfs= addlongorfs.o Sequence.o GeneticCode.o SeqMap.o Array.o
removeoverlaps=removeoverlaps.o Sequence.o SeqMap.o Array.o
scanblastpairs=scanblastpairs.o Sequence.o BlastIndex.o CodonList.o
critica=critica.o Sequence.o GeneticCode.o Statistics.o Scores.o \
        Features.o ScoringMatrix.o DicodonScores.o Regions.o Array.o \
        frameshifts.o 
reportscore=reportscore.o Sequence.o GeneticCode.o Scores.o ScoringMatrix.o \
            DicodonScores.o  
all: lookat motiffind translate extractcoding intergenic empirical-matrix dicodontable addlongorfs removeoverlaps scanblastpairs critica reportscore

lookat: $(lookat)
	$(CC)  -o lookat $(lookat) $(LIBS)

motiffind: $(motiffind)
	$(CC)  -o motiffind $(motiffind) $(LIBS)

translate: $(translate)
	$(CC)  -o translate $(translate) $(LIBS)

extractcoding: $(extractcoding)
	$(CC)  -o extractcoding $(extractcoding) $(LIBS)

intergenic: $(intergenic)
	$(CC)  -o intergenic $(intergenic) $(LIBS)

empirical-matrix: $(empirical-matrix) 
	$(CC)  -o empirical-matrix $(empirical-matrix) $(LIBS) -lm

dicodontable: $(dicodontable) 
	$(CC)  -o dicodontable $(dicodontable) $(LIBS) -lm

addlongorfs: $(addlongorfs) 
	$(CC)  -o addlongorfs $(addlongorfs) $(LIBS) 

removeoverlaps: $(removeoverlaps) 
	$(CC)  -o removeoverlaps $(removeoverlaps) $(LIBS) 

scanblastpairs: $(scanblastpairs) 
	$(CC)  -o scanblastpairs $(scanblastpairs) $(LIBS)

critica: $(critica) 
	$(CC)  -o critica $(critica) $(LIBS) -lm

reportscore: $(reportscore) 
	$(CC)  -o reportscore $(reportscore) $(LIBS) -lm

clean:
	rm -f *.o

debug:
	make clean; make all CC='$(CC) -g'

checker:
	make clean; make all CC='checkergcc -Wall -O3 -g '

efence:
	make clean; make debug LIBS='-lefence'

profile: 
	make clean; make all CC='$(CC) -pg'