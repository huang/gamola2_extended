CC=g++
# CCOPTS=-g -Wall -DDEBUG
CCOPTS=-O2

COGcognitor: main.o cognitor.o enum.o
	$(CC) $(CCOPTS) -o COGcognitor enum.o cognitor.o main.o

main.o: main.cpp cognitor.cpp cognitor.h cognitorglob.h enum.cpp enum.h os.h
	$(CC) $(CCOPTS) -c main.cpp

cognitor.o: cognitor.cpp cognitor.h enum.cpp enum.h
	$(CC) $(CCOPTS) -c cognitor.cpp	

enum.o: enum.cpp enum.h
	$(CC) $(CCOPTS) -c enum.cpp

clean:
	rm *.o COGcognitor
