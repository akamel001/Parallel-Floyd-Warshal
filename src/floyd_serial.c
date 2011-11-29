#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mt19937p.h"
#include <time.h>

#define INF 9999999
#define MIN_RUNS 4
#define MIN_SECS 0.25

const char* usage =
    "path.x -- Parallel all-pairs shortest path on a random graph\n"
    "Flags:\n"
    "  - n -- number of nodes (200)\n"
    "  - p -- probability of including edges (0.05)\n"
    "  - i -- file name where adjacency matrix should be stored (none)\n"
    "  - o -- file name where output matrix should be stored (none)\n";

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

void floyd(int *l, int n){

	for(int k=0; k<n; k++){
			for(int i=0; i<n; i++){
				int ik = i * n + k;
				for (int j=0; j<n; j++){
					int ij = i * n + j;
					int kj = k * n + j;
					if(i == j ) l[ij] = 0;
					if(l[ik]+l[kj]< l[ij])
						l[ij] = l[ik]+l[kj];
				}
		}
	}
}

static inline void infinitize(int n, int* l)
{
    for (int i = 0; i < n*n; ++i)
        if (l[i] == 0)
            l[i] = INF;
}

static inline void deinfinitize(int n, int* l)
{
    for (int i = 0; i < n*n; ++i)
        if (l[i] == INF)
            l[i] = 0;
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

double time_floyd(const int n, const int *l){
	
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
			floyd(lcpy,n);
		finish = clock();
		secs = (finish-start)/(double)(CLOCKS_PER_SEC);
		mflops = 2.0 * num_iterations * n * n * n / 1.0e6; 
		mflop_s = mflops/secs;
		num_iterations *= 2;
		free(lcpy);
	}
	return mflop_s;
}

int main (int argc, char** argv)
{
	int n 	= 200;
	double p = .05;
	const char* ifname = NULL;
	const char* ofname = NULL;
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
	
	 //make graph + output
	 int* l = gen_graph(n,p);
	 if(ifname)
		 write_matrix(ifname, n, l);

	 
	 clock_t start = clock();
	 infinitize(n,l);
	 floyd(l,n); 
	 deinfinitize(n,l);
	 clock_t end = clock();
	 
	 printf("== Serial Floyd Warshall ==\n");
	 printf("n:     %d\n", n);
	 printf("p:     %g\n", p);
	 printf("Time:  %g sec\n", (end-start)/(double)(CLOCKS_PER_SEC));
    if(FLOPS){
	 	mflop_s = time_floyd(n, l);
	 	printf("mflop/s: %lg\n", mflop_s);
	 }
	 printf("Check: %X\n", fletcher16(l, n*n));

	 //output
	 if(ofname)
		 write_matrix(ofname,n,l);

	 free(l);
	 return 0;
}
