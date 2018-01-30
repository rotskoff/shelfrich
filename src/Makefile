CC = g++
DEBUG = -g
GSLINC=#-I/newhome/gmr/gsl-2.3/include
GSL=#-L/newhome/gmr/gsl-2.3/lib
CFLAGS = -Wall $(DEBUG) -O3 -ffast-math -funroll-loops
LFLAGS = -Wall $(DEBUG) $(GSL) -lgsl -lgslcblas -lm

assembly: graph.cc graph.h montecarlo.h montecarlo.cc cargo.cc cargo.h run_assembly.cc
	$(CC) $(CFLAGS) -c graph.cc
	$(CC) $(CFLAGS) -c cargo.cc
	$(CC) $(CFLAGS) -c montecarlo.cc
	$(CC) $(CFLAGS) -c run_assembly.cc
	$(CC) $(CFLAGS) graph.o cargo.o montecarlo.o run_assembly.o -o assemble $(LFLAGS)

umbrella_sampling: graph.cc graph.h montecarlo.h montecarlo.cc cargo.cc cargo.h run_umbrella_sampling.cc
	$(CC) $(CFLAGS) -c graph.cc
	$(CC) $(CFLAGS) -c cargo.cc
	$(CC) $(CFLAGS) -c montecarlo.cc
	$(CC) $(CFLAGS) -c run_umbrella_sampling.cc
	$(CC) $(CFLAGS) graph.o cargo.o montecarlo.o run_umbrella_sampling.o -o umbrella_sampling $(LFLAGS)

testfusion: graph.cc graph.h montecarlo.h montecarlo.cc cargo.cc cargo.h test_fusion.cc
	$(CC) $(CFLAGS) -c graph.cc
	$(CC) $(CFLAGS) -c cargo.cc
	$(CC) $(CFLAGS) -c montecarlo.cc
	$(CC) $(CFLAGS) -c test_fusion.cc
	$(CC) $(CFLAGS) graph.o cargo.o montecarlo.o test_fusion.o -o test_fusion $(LFLAGS)

docs: Doxyfile
	doxygen Doxyfile

clean:
	if [ -f assemble ]; then rm assemble; fi
	if [ -f umbrella_sampling ]; then rm umbrella_sampling; fi
	rm *.o

cleandata:
	if [ -f assemble ]; then rm assemble; fi
	rm lat*.xml
	rm input_parameters.out energy.dat defects.dat

cleandocs:
	if [ -d ../docs ]; then rm -rf docs; fi
