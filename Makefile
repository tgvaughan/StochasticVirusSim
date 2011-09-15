all: hiv_tauleap_hybrid hiv_gillespie

hiv_gillespie: hiv_gillespie.o poissonian.o
	mpic++ -o hiv_gillespie hiv_gillespie.o poissonian.o

hiv_gillespie.o: hiv_gillespie.cc poissonian.h
	mpic++ -c hiv_gillespie.cc

hiv_tauleap_hybrid: hiv_tauleap_hybrid.o poissonian.o
	mpic++ -g -o hiv_tauleap_hybrid hiv_tauleap_hybrid.o poissonian.o

hiv_tauleap_hybrid.o: hiv_tauleap_hybrid.cc poissonian.h
	mpic++ -g -c hiv_tauleap_hybrid.cc

poissonian.o: poissonian.cc poissonian.h
	g++ -c poissonian.cc

clean:
	rm -f hiv_gillespie.o \
		hiv_tauleap_hybrid.o \
		poissonian.o
