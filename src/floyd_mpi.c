#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <algorithm.h>
#include <unistd.h>
#include "mt19937p.h"
#include <time.h>
#include <mpi.h>
#include <assert.h>

char name[MPI_MAX_PROCESSOR_NAME];

#define INF 999999
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

/*
 * Resources:
 * http://www.mcs.anl.gov/~itf/dbpp/text/node35.html
 * http://www.scribd.com/doc/23340655/11/Floyd%E2%80%99s-sequential-algorithm
 */

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

void floyd(int rank, 
			  int *l, 				//data
			  int n, 				//size of matrix
			  int start,			//starting position
			  int block_offset)	
{

	int k,i,j,a;
	int ij,ik;
	//int kj;
	int row_k[n];

	// Maximum path length is N so we iterate N times
	for(k=0; k<n; k++){

		int me = (int) ceil(k/block_offset);

		if(rank == me){
			for(a=0;a<n;a++)
				row_k[a] = l[k*n+a];
		}
		//bcast row 
		MPI_Bcast(&k, 1, MPI_INT, me, MPI_COMM_WORLD);
	
		//bcast column
		MPI_Bcast(row_k, n, MPI_INT, me, MPI_COMM_WORLD);

		for(i=start; i<start+block_offset; i++){
			for (j=0; j<n; j++){
				ij = i * n + j;
				ik = i * n + k;

				if(i==j){
					l[ij] = 0;
				}else{
					//test and set to infinity 
					//infinity redefined since it was defined in math library
					if(l[ij] == 0) l[ij] = INF;
					//TODO write my own min function
					l[ij] = min(l[ij], l[ik]+row_k[j]);
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

void master(int size, int n, int* l){
	MPI_Status status;
	
	//bcast size of matrix
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//bcast the matrix
	MPI_Bcast(l, n*n, MPI_INT, 0, MPI_COMM_WORLD);

	int count = (int) ceil(n/size);
	
	floyd(0,l,n,0,count);

	int recv_buff[n*n];
	for(int p = 1; p < size; p++){
		MPI_Recv(&recv_buff, n*n, MPI_INT, p, 0, MPI_COMM_WORLD, &status);
		for(int i = 0; i < n*n; i++)
			l[i] = max(l[i], recv_buff[i]); //TODO write my own max function
	}

	deinfinitize(n,l);
}

void slave(int rank, int size){
	int n;
	//MPI_Status status;

	//recv the size of matrix
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int lnew[n*n];

	//recv matrix
	MPI_Bcast(&lnew, n*n, MPI_INT, 0, MPI_COMM_WORLD);
	
	int block_offset = (int) ceil(n/size);
	int start = rank * block_offset;
	
	if((n * start) + (n * block_offset) > size) block_offset = n-start;

	floyd(rank, lnew, n, start, block_offset);
	
	//send back my data to master processor 
	MPI_Send(lnew, n*n, MPI_INT, 0,0, MPI_COMM_WORLD);
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
	int processors, my_rank, size;
    int* my_partition;
    int i, j, k, partition_size;
    int* kth_row;
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
	int* graph = gen_graph(n,p);
	if(ifname)
		write_matrix(ifname, n, graph);
    size = n;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processors);
    

    /*----------
     // Broadcast data size.
     //--------*/
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    /*----------
     // Determine data partitions.
     //--------*/
    partition_size = size / processors;
    if (size % processors)
        partition_size++;
    
    my_partition = (int*) malloc(size * partition_size * sizeof(int));
    kth_row      = (int*) malloc(size * 1 * sizeof(int));
    
    assert(my_partition != NULL);
    
    /*----------
     // Distribute the data.
     //--------*/
    MPI_Scatter(graph, partition_size * size, MPI_INT, my_partition,
                partition_size * size, MPI_INT, 0, MPI_COMM_WORLD);
    
    /* If processors does not divide evenly into size, the Scatter
     will send random memory beyond the end of graph, but that data
     should not actually get used, so it probably won't crash. */
    
    /*----------
     // Do the calculation.
     //--------*/
    for (k = 0; k < size; k++) {
        
        /*----------
         // Broadcast kth row.
         //--------*/
        if (my_rank == (k / partition_size))
            for (i = 0; i < size; i++)
                kth_row[i] = my_partition[(k%partition_size)*size+i];
        MPI_Bcast(kth_row, size, MPI_INT, (k / partition_size), MPI_COMM_WORLD);
        
        /*----------
         // Update my rows.
         //--------*/
        for (i = 0; (i < partition_size) && (i < size); i++) {
            if (my_partition[i*size+k] < 1)
                continue;
            for (j = 0; j < size; j++) {
                if (kth_row[j] < 1)
                    continue;
                if (my_partition[i*size+j] < 0)
                    my_partition[i*size+j] = my_partition[i*size+k] + kth_row[j];
                else
                    my_partition[i*size+j] = min(my_partition[i*size+j],
                                                 my_partition[i*size+k] + kth_row[j]);
            }
        }
    }
    
    /*----------
     // Collect the data.
     //--------*/
    MPI_Gather(my_partition, partition_size * size, MPI_INT, graph,
               partition_size * size, MPI_INT, 0, MPI_COMM_WORLD);

//	MPI_Get_processor_name(name, &len);//TODO is this needed?
//
	if(my_rank == 0){
//		clock_t start = clock();
//		master(size, n, l);
//		clock_t end = clock();
//		
		printf("== MPI Floyd-Warshall %d processors\n", processors);
		printf("n:     %d\n", n);
		printf("p:     %g\n", p);
		//printf("Time:  %g us\n", (double) (end-start)/1000);
		printf("Check: %X\n", fletcher16(graph, n*n));
	
		//output
		if(ofname)
			write_matrix(ofname,n,graph);
	} 
//	else {
//		slave(rank, size);
//	}

	//infinitize(n,l);
	//floyd(l,n); 
	//deinfinitize(n,l);


	MPI_Finalize();
	free(graph);
	return 0;
}






