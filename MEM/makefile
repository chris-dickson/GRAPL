CC=gcc

#CFLAGS= -g -c -fnested-functions -O2
#LDFLAGS= -g -fnested-functions -O2

#Dantzig
CFLAGS=-O2 -g -c
LDFLAGS=-O2 -g

BUILDDIR=./build/

all :  route frac_round_MRD prune_yal prune_ibm route_MT

route : graph.o list.o set.o helper.o graphtype.o linesearch.o parseargs.o rounding.o routing.o steiner.o mst.o shortest.o binheap.o main.c globals.o 
	$(CC) $(LDFLAGS) -lm -lpthread graph.o list.o set.o helper.o graphtype.o linesearch.o parseargs.o rounding.o routing.o steiner.o mst.o shortest.o binheap.o globals.o main.c -o route
	cp route $(BUILDDIR)
	
route_MT : graph.o list.o set.o helper.o graphtype.o linesearch.o parseargs.o rounding.o routing.o steiner.MT.o mst.o shortest.MT.o binheap.o main.MT.c globals.MT.o
	$(CC) $(LDFLAGS) -lm -lpthread graph.o list.o set.o helper.o graphtype.o linesearch.o parseargs.o rounding.o routing.o steiner.MT.o mst.o shortest.MT.o binheap.o globals.MT.o main.MT.c -o route_MT
	cp route_MT $(BUILDDIR)

route_ECMM : graph.o list.o set.o helper.o graphtype.o linesearch.o parseargs.o rounding.o routing.o steiner.o mst.o shortest.o binheap.o main.ecmm.c globals.o 
	$(CC) $(LDFLAGS) -lm graph.o list.o set.o helper.o graphtype.o linesearch.o parseargs.o rounding.o routing.o steiner.o mst.o shortest.o binheap.o globals.o main.ecmm.c -o route_ECMM
	cp route_ECMM $(BUILDDIR)

	
frac_round_MRD: graph.o list.o set.o helper.o rounding.o frac_round_MRD.c frac_round_aux.c
	$(CC) $(LDFLAGS) -lm graph.o list.o set.o helper.o rounding.o frac_round_MRD.c -o frac_round_MRD
	cp frac_round_MRD $(BUILDDIR)
	
prune_yal: list.o set.o prune_yal.c
	$(CC) $(LDFLAGS) list.o set.o prune_yal.c -o prune_yal
	cp prune_yal $(BUILDDIR)

prune_ibm: list.o set.o prune_ibm.c
	$(CC) $(LDFLAGS) list.o set.o prune_ibm.c -o prune_ibm
	cp prune_ibm $(BUILDDIR)
	
steiner.o : steiner.h steiner.c graph.o list.o mst.o shortest.o routing.o
	$(CC) $(CFLAGS) steiner.c -o steiner.o
	
steiner.MT.o : steiner.MT.h steiner.MT.c graph.o list.o mst.o shortest.o routing.h
	$(CC) $(CFLAGS) steiner.MT.c -o steiner.MT.o

rounding.o : rounding.h rounding.c
	$(CC) $(CFLAGS) rounding.c -o rounding.o

linesearch.o : linesearch.h linesearch.c 
	$(CC) $(CFLAGS) linesearch.c -o linesearch.o
	
shortest.o : shortest.h shortest.c binheap.o graph.o list.o set.o helper.o globals.o
	$(CC) $(CFLAGS) shortest.c -o shortest.o
	
shortest.MT.o : shortest.MT.h shortest.MT.c binheap.o graph.o list.o set.o helper.o globals.MT.o
	$(CC) $(CFLAGS) shortest.MT.c -o shortest.MT.o

mst.o : mst.h mst.c graph.o list.o
	$(CC) $(CFLAGS) mst.c -o mst.o

graphtype.o : graphtype.h graphtype.c graph.o list.o set.o helper.o
	$(CC) $(CFLAGS) graphtype.c -o graphtype.o

routing.o : routing.h routing.c helper.o set.o graph.o
	$(CC) $(CFLAGS) routing.c -o routing.o

binheap.o : binheap.h binheap.c
	$(CC) $(CFLAGS) binheap.c -o binheap.o

parseargs.o : parseargs.h parseargs.c
	$(CC) $(CFLAGS) parseargs.c -o parseargs.o

helper.o : helper.h helper.c graph.o list.o set.o
	$(CC) $(CFLAGS) helper.c -o helper.o
	
globals.o : globals.c globals.h binheap.o helper.h
	$(CC) $(CFLAGS) globals.c -o globals.o
	
globals.MT.o : globals.MT.c globals.MT.h binheap.o helper.h
	$(CC) $(CFLAGS) globals.MT.c -o globals.MT.o

graph.o : graph.c graph.h graphalg.h list.o set.o
	$(CC) $(CFLAGS) graph.c -o graph.o

set.o : set.h set.c list.h
	$(CC) $(CFLAGS) set.c -o set.o

list.o : list.h list.c
	$(CC) $(CFLAGS) list.c -o list.o

clean : 
	rm -f *.o
	rm -f $(BUILDDIR)route $(BUILDDIR)frac_round_MRD $(BUILDDIR)route_MT
	rm -f $(BUILDDIR)AllTrees.txt $(BUILDDIR)Stats.txt
	rm -r -f $(BUILDDIR)output/
	rm -f $(BUILDDIR)prune_*
	rm -f route frac_round_MRD route_MT prune_yal prune_ibm route_ECMM

archive :
	tar cvzf archive.tgz $(BUILDDIR)AllTrees.txt $(BUILDDIR)Stats.txt $(BUILDDIR)output/*
