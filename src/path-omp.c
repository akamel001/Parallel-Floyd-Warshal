#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <time.h>
#include "mt19937p.h"

#define MIN_RUNS 4
#define MIN_SECS 0.25

/*
 * Description: Algorithm uses repeated squaring process to solve the matrix
 */
int square(int n,               // Number of nodes
           int* restrict l,     // Partial distance at step s
           int* restrict lnew)  // Partial distance at step s+1
{
    int done = 1;
    #pragma omp parallel for shared(l, lnew) reduction(&& : done)
    for (int j = 0; j < n; ++j) {
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

void shortest_paths(int n, int* restrict l)
{
    infinitize(n, l);
    for (int i = 0; i < n*n; i += n+1)
        l[i] = 0;

    // Repeated squaring until nothing changes
    int* restrict lnew = (int*) calloc(n*n, sizeof(int));
    memcpy(lnew, l, n*n * sizeof(int));
	 for (int done = 0; !done; ) {
        done = square(n, l, lnew);
		 memcpy(l, lnew, n*n * sizeof(int));
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

double time_path(const int n, const int *l){
	
	clock_t start, finish;
	double mflops, mflop_s;
	double secs = -1.0;
	int num_iterations = MIN_RUNS;
	int i;
    
	while(secs < MIN_SECS){
		int *lcpy = (int*) calloc(n*n, sizeof(int));
		memcpy(lcpy,l, n*n*sizeof(int));
		start = clock();
		for(i = 0; i < num_iterations; ++i)
			shortest_paths(n, lcpy);
		finish = clock();
		secs = (finish-start)/(double)(CLOCKS_PER_SEC);
		mflops = 2.0 * num_iterations * n * n * n / 1.0e6; 
		mflop_s = mflops/secs;
		num_iterations *= 2;
		free(lcpy);
	}
	return mflop_s;
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
    double mflop_s;
    int FLOPS = 0;
    
    //args
	extern char* optarg;
	const char* optstring = "hn:d:p:o:i:f:";
	int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch (c) {
            case 'h':
                fprintf(stderr, "%s", usage);
                return -1;
            case 'n': n = atoi(optarg);			break;
            case 'p': p = atof(optarg); 		break;
            case 'o': ofname = optarg;  		break;
            case 'f': FLOPS = atoi(optarg); 	break;
            case 'i': ifname = optarg; 		 	break;
        }
    }

    // Graph generation + output
    int* l = gen_graph(n, p);
    if (ifname)
        write_matrix(ifname,  n, l);

    // Time the shortest paths code
    double t0 = omp_get_wtime();
    shortest_paths(n, l);
    double t1 = omp_get_wtime();

    printf("== OpenMP with %d threads\n", omp_get_max_threads());
    printf("n:     %d\n", n);
    printf("p:     %g\n", p);
    printf("Time:  %g\n", t1-t0);
    if(FLOPS){
	 	mflop_s = time_path(n, l);
	 	printf("mflop/s: %lg\n", mflop_s);
    }
    printf("Check: %X\n", fletcher16(l, n*n));

    // Generate output file
    if (ofname)
        write_matrix(ofname, n, l);

    // Clean up
    free(l);
    return 0;
}
