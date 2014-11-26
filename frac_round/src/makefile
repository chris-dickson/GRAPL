all : frac_round_MRD frac_round_TARGET
frac_round_TARGET : graph.c graph.h graphalg.h list.c list.h main.mrd.target.c rounding.c set.c set.h
	gcc -g main.mrd.target.c -o frac_round_TARGET
frac_round_MRD : graph.c graph.h graphalg.h list.c list.h main.mrd.c rounding.c set.c set.h
	gcc -g main.mrd.c -o frac_round_MRD
frac_round_OVR : graph.c graph.h graphalg.h list.c list.h main.ovr.c rounding.c set.c set.h
	gcc -g main.ovr.c -o frac_round_OVR
clean : 
	rm frac_round_MRD frac_round_TARGET
