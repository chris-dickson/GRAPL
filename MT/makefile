CC=gcc
CFLAGS= -g -c -O2
LDFLAGS= -g -O2
BUILDDIR=./build/
MEMDIR=../MEM/

all :  route_MT frac_round_MRD prune_yal prune_ibm

route_MT : graph.o list.o set.o helper.o graphtype.o linesearch.o parseargs.o rounding.o routing.o steiner.o mst.o shortest.o binheap.o globals.o main.c
	$(CC) $(LDFLAGS) -lpthread -lm graph.o list.o set.o helper.o graphtype.o linesearch.o parseargs.o rounding.o routing.o steiner.o mst.o shortest.o binheap.o globals.o main.c -o route_MT
	mv route_MT $(BUILDDIR)


frac_round_MRD: $(MEMDIR)graph.o $(MEMDIR)list.o $(MEMDIR)set.o $(MEMDIR)helper.o $(MEMDIR)rounding.o $(MEMDIR)frac_round_MRD.c $(MEMDIR)frac_round_aux.c
	$(CC) $(LDFLAGS) -lm $(MEMDIR)graph.o $(MEMDIR)list.o $(MEMDIR)set.o $(MEMDIR)helper.o $(MEMDIR)rounding.o $(MEMDIR)frac_round_MRD.c -o frac_round_MRD
	mv frac_round_MRD $(BUILDDIR)
	
prune_yal: $(MEMDIR)list.o $(MEMDIR)set.o $(MEMDIR)prune_yal.c
	$(CC) $(LDFLAGS) $(MEMDIR)list.o $(MEMDIR)set.o $(MEMDIR)prune_yal.c -o prune_yal
	mv prune_yal $(BUILDDIR)

prune_ibm: $(MEMDIR)list.o $(MEMDIR)set.o $(MEMDIR)prune_ibm.c
	$(CC) $(LDFLAGS) $(MEMDIR)list.o $(MEMDIR)set.o $(MEMDIR)prune_ibm.c -o prune_ibm
	mv prune_ibm $(BUILDDIR)
	
steiner.o : steiner.h steiner.c graph.o list.o mst.o shortest.o routing.o
	$(CC) $(CFLAGS) steiner.c -o steiner.o

rounding.o : rounding.h rounding.c
	$(CC) $(CFLAGS) rounding.c -o rounding.o

linesearch.o : linesearch.h linesearch.c 
	$(CC) $(CFLAGS) linesearch.c -o linesearch.o
	
shortest.o : shortest.h shortest.c binheap.o graph.o list.o set.o helper.o globals.o
	$(CC) $(CFLAGS) shortest.c -o shortest.o

mst.o : mst.h mst.c graph.o list.o
	$(CC) $(CFLAGS) mst.c -o mst.o

graphtype.o : graphtype.h graphtype.c graph.o list.o set.o helper.o
	$(CC) $(CFLAGS) graphtype.c -o graphtype.o

routing.o : routing.h routing.c helper.o set.o
	$(CC) $(CFLAGS) routing.c -o routing.o

binheap.o : binheap.h binheap.c
	$(CC) $(CFLAGS) binheap.c -o binheap.o

parseargs.o : parseargs.h parseargs.c
	$(CC) $(CFLAGS) parseargs.c -o parseargs.o

helper.o : helper.h helper.c graph.o list.o set.o
	$(CC) $(CFLAGS) helper.c -o helper.o
	
globals.o : globals.c globals.h binheap.o graph.o graphtype.o routing.o
	$(CC) $(CFLAGS) globals.c -o globals.o

graph.o : graph.c graph.h graphalg.h list.o set.o
	$(CC) $(CFLAGS) graph.c -o graph.o

set.o : set.h set.c list.h
	$(CC) $(CFLAGS) set.c -o set.o

list.o : list.h list.c
	$(CC) $(CFLAGS) list.c -o list.o

clean : 
	rm -f *.o
	rm -f $(BUILDDIR)route_MT
	rm -f $(BUILDDIR)frac_round_MRD
	rm -f $(BUILDDIR)prune_yal $(BUILDDIR)prune_ibm
	rm -f $(BUILDDIR)AllTrees.txt $(BUILDDIR)Stats.txt
	rm -r -f $(BUILDDIR)output/

archive :
	tar cvzf archive.tgz $(BUILDDIR)AllTrees.txt $(BUILDDIR)Stats.txt $(BUILDDIR)output/*
