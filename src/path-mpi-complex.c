#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mpi.h"
#include "mt19937p.h"

/*
 * Testing and Sanity Check purpose while debugging
 */
void write_matrix2(const char* fname, int n, int nloc, int* a)
{
	FILE* fp = fopen(fname, "w+");
	if (fp == NULL) {
		fprintf(stderr, "Could not open output file: %s\n", fname);
		exit(-1);
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < nloc; ++j){ 
			fprintf(fp, "%d ", a[j*n+i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

/* --
 * Exchange ghost cell data with neighboring processors
 */
void exchange(int rank, int size, int count, int* pass_buff, int* tmp_buff)
{
	if (size == 1)
		return;

	if (rank == size-1) {
		MPI_Sendrecv(pass_buff, count, MPI_INT, rank-1, 0, tmp_buff,
				count, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
	} else if (rank == 0) {
		MPI_Sendrecv(pass_buff, count, MPI_INT, size-1, 0, tmp_buff,
				count, MPI_INT, rank+1, 0, MPI_COMM_WORLD, NULL);
	} else {
		MPI_Sendrecv(pass_buff, count, MPI_INT, rank-1, 0, tmp_buff,
				count, MPI_INT, rank+1, 0, MPI_COMM_WORLD, NULL);
	}
	memcpy(pass_buff, tmp_buff, count*sizeof(int));
}

int square(int n,               // Number of nodes
		int nloc,
		int* restrict lloc,     // Partial distance at step s
		int* restrict lnew,  // Partial distance at step s+1
		int* restrict pass_buff,
		int* tmp_buff,
		int block_size,
		int rank,
		int size)
{
	int done = 1;
	for (int phase = 0; phase < size; ++phase) {
		if (phase > 0) {
			exchange(rank, size, block_size*n, pass_buff, tmp_buff);
		}
		int br = (rank + phase) % size;
		int bi = br * block_size;
		int kMax = (br == size-1) ? n - bi : block_size;
		for (int j = 0; j < nloc; ++j) {
			for (int i = 0; i < n; ++i) {
				int lij = lnew[j*n+i];
				for (int k = 0; k < kMax; ++k) {
					int lik = pass_buff[k*n+i];
					int lkj = lloc[j*n+k+bi];
					if (lik + lkj < lij) {
						lij = lik+lkj;
						done = 0;
					}
				}
				lnew[j*n+i] = lij;
			}
		}
	}
	memcpy(lloc, lnew, n*nloc * sizeof(int));
	MPI_Allreduce(&done, &done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
	return done;
}

static inline void infinitize(int n, int nloc, int* l)
{
	for (int i = 0; i < n*nloc; ++i)
		if (l[i] == 0)
			l[i] = n+1;
}

static inline void deinfinitize(int n, int nloc, int* l)
{
	for (int i = 0; i < n*nloc; ++i)
		if (l[i] == n+1)
			l[i] = 0;
}

void shortest_paths(int n,
		int nloc,
		int* restrict lloc,
		int* restrict lnew,
		int* restrict pass_buff,
		int* restrict tmp_buff,
		int block_size,
		int rank,
		int size)
{
	infinitize(n, nloc, lloc);
	for (int i = rank*block_size; i < n*nloc; i += n+1)
		lloc[i] = 0;

	memcpy(lnew, lloc, n*nloc * sizeof(int));
	memset(pass_buff, n+1, n*block_size*sizeof(int));
	memset(tmp_buff, n+1, n*block_size*sizeof(int));
	int done = 0;
	for (done = 0; !done; ) {
		memcpy(pass_buff, lloc, n*nloc * sizeof(int));
		done = square(n, nloc, lloc, lnew,
				pass_buff, tmp_buff, block_size, rank, size);
	}
	deinfinitize(n, nloc, lloc);

}

int* gen_graph(int n, double p)
{
	int* l = calloc(n*n, sizeof(int));
	struct mt19937p state;
	sgenrand(10302011UL, &state);
	for (int j = 0; j < n; ++j) {
		for (int i = 0; i < n; ++i)
			l[j*n+i] = (genrand(&state) < p);
		l[j*n+j] = 0;
	}
	return l;
}

int fletcher16(int* data, int count)
{
	int sum1 = 0;
	int sum2 = 0;
	for(int index = 0; index < count; ++index) {
		sum1 = (sum1 + data[index]) % 255;
		sum2 = (sum2 + sum1) % 255;
	}
	return (sum2 << 8) | sum1;
}

void write_matrix(const char* fname, int n, int* a)
{
	FILE* fp = fopen(fname, "w+");
	if (fp == NULL) {
		fprintf(stderr, "Could not open output file: %s\n", fname);
		exit(-1);
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) 
			fprintf(fp, "%d ", a[j*n+i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

const char* usage =
"path.x -- Parallel all-pairs shortest path on a random graph\n"
"Flags:\n"
"  - n -- number of nodes (200)\n"
"  - p -- probability of including edges (0.05)\n"
"  - i -- file name where adjacency matrix should be stored (none)\n"
"  - o -- file name where output matrix should be stored (none)\n";

int main(int argc, char** argv)
{
	int n    = 200;            // Number of nodes
	double p = 0.05;           // Edge probability
	const char* ifname = NULL; // Adjacency matrix file name
	const char* ofname = NULL; // Distance matrix file name
	int rank, size;
	double t0, t1;

	// Option processing
	extern char* optarg;
	const char* optstring = "hn:d:p:o:i:";
	int c;
	while ((c = getopt(argc, argv, optstring)) != -1) {
		switch (c) {
			case 'h':
				fprintf(stderr, "%s", usage);
				return -1;
			case 'n': n = atoi(optarg); break;
			case 'p': p = atof(optarg); break;
			case 'o': ofname = optarg;  break;
			case 'i': ifname = optarg;  break;
		}
	}

	// Initialize MPI and get rank and size
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int block_size, joffset, nloc;

	// Split internal index range among processors
	block_size    = ceil((double)n / (double)size);
	joffset       = rank * block_size;
	nloc          = (rank == size-1) ? (n-joffset) : block_size;

	// Graph generation + output
	int* l = gen_graph(n, p);
	if (rank==0 && ifname)
		write_matrix(ifname,  n, l);

	// Repeated squaring until nothing changes
	int* restrict lnew      = (int*) calloc(n*nloc, sizeof(int));
	int* restrict pass_buff = (int*) calloc(n*block_size, sizeof(int));
	int* restrict tmp_buff  = (int*) calloc(n*block_size, sizeof(int));

	// Allocate and initialize local array
	int* lloc = (int*) calloc( nloc*n, sizeof(int) );
	memcpy(lloc, &l[joffset*n], n*nloc*sizeof(int));

	// Time the shortest paths code
	t0 = MPI_Wtime();
	shortest_paths(n, nloc, lloc, lnew, pass_buff, tmp_buff, block_size, rank, size);
	t1 = MPI_Wtime();

	int displs[size], recv_counts[size];	
	for (int r = 0; r < size; ++r) {
		displs[r] = r * block_size * n;
		recv_counts[r] =  block_size * n;
	}

	recv_counts[size-1] = n*n - (size-1)*block_size*n;
	MPI_Allgatherv(lloc, n*nloc, MPI_INT, l,
			recv_counts, displs, MPI_INT, MPI_COMM_WORLD);
	
	if(rank == 0){	
		printf("== MPI Complex with %d processes\n", size);
		printf("n:     %d\n", n);
		printf("p:     %g\n", p);
		printf("Time:  %g\n", t1-t0);
		printf("Check: %X\n", fletcher16(l, n*n));

		// Generate output file
		if (ofname)
			write_matrix(ofname, n, l);
	}

	// Clean up
	free(lnew);
	free(pass_buff);
	free(tmp_buff);
	free(lloc);
	free(l);
	MPI_Finalize();
	return 0;
}
