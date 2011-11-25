#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mt19937p.h"
#include <time.h>

#define INF 999999

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

void floyds(int *l, int n){

	int k,i,j;
	int ij,ik,kj;

	// Maximum path length is N so we iterate N times
	for(k=0; k<n; k++){
		for(i=0; i<n; i++){
			ik = i * n + k;
			for (j=0; j<n; j++){
				ij = i * n + j;
				kj = k * n + j;
				if(i==j){
					l[ij] = 0;
				}else{
					//test and set to infinity 
					//infinity redefined since it was defined in math library
					if(l[ij] == 0) l[ij] = INF;
					
					// If our data is smaller, replace it 
					// and set the output to be the current path length
					if(l[ik]+l[kj]< l[ij]){
						l[ij] = l[ik]+l[kj];
					}
				}
			}
		}
	}
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

int main (int argc, char** argv)
{
	int n 	= 200;
	double p = .05;
	const char* ifname = NULL;
	const char* ofname = NULL;

	//args
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
	
	 //make graph + output
	 int* l = gen_graph(n,p);
	 if(ifname)
		 write_matrix(ifname, n, l);

	 clock_t start = clock();
	 //infinitize(n,l);
	 floyds(l,n); 
	 deinfinitize(n,l);
	 clock_t end = clock();
	
    printf("== Serial FloydWarshall\n");
    printf("n:     %d\n", n);
    printf("p:     %g\n", p);
    printf("Time:  %g\n", (double) (end-start));
    printf("Check: %X\n", fletcher16(l, n*n));
	
	//output
	if(ofname)
		write_matrix(ofname,n,l);

	free(l);
	return 0;
}






