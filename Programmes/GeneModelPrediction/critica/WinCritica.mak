CC = $(CC) /G6 /W3 /DMICROSOFT

lookat=lookat.obj Sequence.obj
motiffind=motiffind.obj Sequence.obj
translate=translate.obj Sequence.obj GeneticCode.obj
extractcoding=extract-coding.obj Sequence.obj GeneticCode.obj 
intergenic=intergenic.c Sequence.obj SeqMap.obj
empiricalmatrix=empiricalmatrix.obj Sequence.obj GeneticCode.obj
dicodontable=dicodontable.obj Sequence.obj SeqMap.obj
addlongorfs= addlongorfs.obj Sequence.obj GeneticCode.obj SeqMap.obj Array.obj
removeoverlaps=removeoverlaps.obj Sequence.obj SeqMap.obj Array.obj
scanblastpairs=scanblastpairs.obj Sequence.obj BlastIndex.obj CodonList.obj
critica=critica.obj Sequence.obj GeneticCode.obj Statistics.obj Scores.obj \
        Features.obj ScoringMatrix.obj DicodonScores.obj Regions.obj Array.obj \
        frameshifts.obj 
reportscore=reportscore.obj Sequence.obj GeneticCode.obj Scores.obj ScoringMatrix.obj \
            DicodonScores.obj  
all: lookat.exe motiffind.exe translate.exe extractcoding.exe intergenic.exe empiricalmatrix.exe dicodontable.exe addlongorfs.exe removeoverlaps.exe scanblastpairs.exe critica.exe reportscore.exe

lookat.exe: $(lookat)
	$(CC)  -o lookat.exe $(lookat) 


motiffind.exe: $(motiffind)
	$(CC)  -o motiffind.exe $(motiffind) 

translate.exe: $(translate)
	$(CC)  -o translate.exe $(translate) 

extractcoding.exe: $(extractcoding)
	$(CC)  -o extractcoding.exe $(extractcoding) 

intergenic.exe: $(intergenic)
	$(CC)  -o intergenic $(intergenic) 

empiricalmatrix.exe: $(empiricalmatrix) 
	$(CC)  -o empiricalmatrix.exe $(empiricalmatrix)  

dicodontable.exe: $(dicodontable) 
	$(CC)  -o dicodontable.exe $(dicodontable)  

addlongorfs.exe: $(addlongorfs) 
	$(CC)  -o addlongorfs.exe $(addlongorfs)  

removeoverlaps.exe: $(removeoverlaps) 
	$(CC)  -o removeoverlaps.exe $(removeoverlaps)  

scanblastpairs.exe: $(scanblastpairs) 
	$(CC)  -o scanblastpairs.exe $(scanblastpairs) 

critica.exe: $(critica) 
	$(CC)  -o critica.exe $(critica)  

reportscore.exe: $(reportscore) 
	$(CC)  -o reportscore.exe $(reportscore)  


