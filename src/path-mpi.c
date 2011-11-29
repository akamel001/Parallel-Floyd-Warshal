#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mpi.h"
#include "mt19937p.h"


/*
 * Description: Algorithm uses repeated squaring process to solve the matrix
 */
int square(int n,          // Number of nodes
		int* restrict l,     // Partial distance at step s
		int* restrict lnew,  // Partial distance at step s+1
		int rank,
		int size)
{
	int done = 1;

	int block_size = n / size; 
	int begin_row = rank * block_size;
	int send_count = block_size * n;
	int newMax; 
	
	if (rank == size-1) {
		send_count = n*n - begin_row*n;
		newMax = n;
	} else 
		newMax = begin_row+block_size;

	int displs[size];
	int recv_counts[size];
	
	for (int j = begin_row; j < newMax; ++j) {
		for (int i = 0; i < n; ++i) {
			int lij = lnew[j*n+i];
			for (int k = 0; k < n; ++k) {
				int lik = l[k*n+i];
				int lkj = l[j*n+k];
				if (lik + lkj < lij) {
					lij = lik+lkj;
					done = 0;
				}
			}
			lnew[j*n+i] = lij;
		}
	}
	
	for (int r = 0; r < size; ++r) {
		displs[r] = r * block_size * n;
		recv_counts[r] =  block_size * n;
	}
	
	recv_counts[size-1] = n*n - (size-1)*block_size*n;

	MPI_Allgatherv(&lnew[begin_row*n], send_count, MPI_INT, l,
			recv_counts, displs, MPI_INT, MPI_COMM_WORLD);
	
	MPI_Allreduce(&done, &done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
	
	
	return done;
}

static inline void infinitize(int n, int* l)
{
	for (int i = 0; i < n*n; ++i)
		if (l[i] == 0)
			l[i] = n+1;
}

static inline void deinfinitize(int n, int* l)
{
	for (int i = 0; i < n*n; ++i)
		if (l[i] == n+1)
			l[i] = 0;
}

void shortest_paths(int n, int* restrict l, int rank, int size)
{
	infinitize(n, l);
	for (int i = 0; i < n*n; i += n+1)
		l[i] = 0;

	// Repeated squaring until nothing changes
	int* restrict lnew = (int*) calloc(n*n, sizeof(int));
	memcpy(lnew, l, n*n * sizeof(int));
	int done = 0;
	for (done = 0; !done; ) {
		done = square(n, l, lnew, rank, size);
	}
	free(lnew);
	deinfinitize(n, l);
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

/* http://en.wikipedia.org/wiki/Fletcher's_checksum} */

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

    // Graph generation + output
    int* l = gen_graph(n, p);
    if (ifname)
        write_matrix(ifname,  n, l);

    // Time the shortest paths code
    t0 = MPI_Wtime();
    shortest_paths(n, l, rank, size);
    t1 = MPI_Wtime();

    if (rank == 0) {
        printf("== MPI with %d processes ==\n", size);
        printf("n:     %d\n", n);
        printf("p:     %g\n", p);
        printf("Time:  %g\n", t1-t0);
        printf("Check: %X\n", fletcher16(l, n*n));
    }

    // Generate output file
    if (ofname)
        write_matrix(ofname, n, l);

    // Clean up
    free(l);
    
    MPI_Finalize();
    return 0;
}
