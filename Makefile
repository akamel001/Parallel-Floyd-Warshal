CC=gcc
MPICC=mpicc
CFLAGS=-Wall -falign-functions -falign-loops -funroll-all-loops -std=gnu99 -O3
DEBUG=no 

ifeq ($(DEBUG),yes)
	CFLAGS=CFLAGS + -g
endif

INCPATH=inc
SRCPATH=src
OBJPATH=obj
BINPATH=bin

INCLUDES=-I ./$(INCPATH)

.PHONY: exe clean realclean


#==========
exe: $(BINPATH)/path-omp \
	  $(BINPATH)/path-mpi \
	  $(BINPATH)/path-mpi-complex \
	  $(BINPATH)/floyd_serial \
	  $(BINPATH)/floyd_omp 

$(BINPATH)/path-omp: $(OBJPATH)/path-omp.o $(OBJPATH)/mt19937p.o
	$(CC) -fopenmp $(CFLAGS) -o $@ $^
 
$(BINPATH)/floyd_omp: $(OBJPATH)/floyd_omp.o $(OBJPATH)/mt19937p.o
	$(CC) -fopenmp $(CFLAGS) -o $@ $^ 

$(BINPATH)/floyd_serial: $(OBJPATH)/floyd_serial.o $(OBJPATH)/mt19937p.o
	$(CC) $(CFLAGS) -o $@ $^ 

$(BINPATH)/path-mpi: $(OBJPATH)/path-mpi.o $(OBJPATH)/mt19937p.o
	$(MPICC) $(CFLAGS) -o $@ $^ 

$(BINPATH)/path-mpi-complex: $(OBJPATH)/path-mpi-complex.o $(OBJPATH)/mt19937p.o
	$(MPICC) $(CFLAGS) -o $@ $^ 

$(OBJPATH)/path-omp.o: $(SRCPATH)/path-omp.c 
	$(CC) -c -fopenmp $(CFLAGS) $(INCLUDES) $< -o $@

$(OBJPATH)/floyd_omp.o: $(SRCPATH)/floyd_omp.c 
	$(CC) -c -fopenmp $(CFLAGS) $(INCLUDES) $< -o $@

$(OBJPATH)/floyd_serial.o: $(SRCPATH)/floyd_serial.c 
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(OBJPATH)/path-mpi.o: $(SRCPATH)/path-mpi.c 
	$(MPICC) -c  $(CFLAGS) $(INCLUDES) $< -o $@

$(OBJPATH)/path-mpi-complex.o: $(SRCPATH)/path-mpi-complex.c 
	$(MPICC) -c  $(CFLAGS) $(INCLUDES) $< -o $@

$(OBJPATH)/mt19937p.o: $(SRCPATH)/mt19937p.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(OBJPATH)/%.o: $(SRCPATH)/%.c
	$(CC) $(INCLUDES) -c $(CFLAGS) $< $@

# === Cleanup and tarball

clean:
	rm -f $(OBJPATH)/*.o

realclean: clean
	rm -f $(BINPATH)/* 

tar:
	(tar -cvzf final.tgz ../final/)
