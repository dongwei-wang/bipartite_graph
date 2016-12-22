#include <iostream>
#include <utility>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <time.h>
//#include "graphfilter.h"

//#include <cutil.h>

//#include <graphfilter_kernel.cu>
#define N 400
#define DEBUG 0

typedef struct _Graph_node_A {
	int adj[N];
} Graph_node_A;

using namespace std;

extern "C"

#ifndef _MATRIXMUL_H_
#define _MATRIXMUL_H_
__global__ void naiveFilterKernel(struct _Graph_node_A *PA, int *relation_count, int *v_start, int *max_id_A)
{
	__shared__ int v_start_adj[N];

	int tid = threadIdx.x;
	int count=0;

	if(tid < N)
		v_start_adj[tid] = PA[(*v_start)-1].adj[tid];
	__syncthreads();

	if(tid == ((*v_start)-1) || tid >= *max_id_A) return;

	int i, j;

	for(i=0;(i<N) && (PA[tid].adj[i] != 0);i++) {
		for(j=0;j<N && v_start_adj[j]!= 0;j++) {
			if(PA[tid].adj[i] == v_start_adj[j]){
				count++;
				break;
			}
		}
	}
	relation_count[tid] = count;
}
#endif

#ifndef _MATRIXMUL_F_
#define _MATRIXMUL_F_
__global__ void fullFilterKernel(struct _Graph_node_A *PA, int *relation_count, int *index_vertice, int *k, int *max_id_A)
{
	__shared__ int index_ref[N];

	int by = blockIdx.y;
	int tx = threadIdx.x;
	int count = 0;

	if((by >= *k) || (tx >= *max_id_A)) return;

	// load data into shared memory
	if(tx < N)
		index_ref[tx] = PA[index_vertice[by]-1].adj[tx];
	__syncthreads();

	for(int i=0;(i<N) && (PA[tx].adj[i] != 0);i++)
		for(int j=0;(j<N) && (index_ref[j]!= 0);j++) {
			if((PA[tx].adj[i] == index_ref[j]) && ((index_vertice[by]-1) != tx)) {
				count++;
				break;
			}
		}

	relation_count[(by+1) * (*max_id_A) + tx] = count;
}
#endif

void filterOnDevice(struct _Graph_node_A *PA, int *relation_count, int *index_vertice, int v_start, int k, int max_id_A);
void filterOnHost(struct _Graph_node_A *PA, int *relation_count_h, int *index_vertice_h, int v_start, int k, int max_id_A);

// function to sort the relations to the start vertice in a descending order
void selSort(int s[], int index[], int length);

//function to compare the result between the host side and device side
bool compare(vector<int> final_resut_h, vector<int> final_result_d);

//static const char* DAT_FILE_NAME = "graph.dat";// file to store the graph struct

int main(int arcg, char** argv)
{
	char tracefilename[30];
	int A, B; //fot the nodes in party A and party B
	int max_id_A;
	int max_id_B;
	int max_weight;
	int i=0, j=0;
	//int m=0, n=0;
	FILE *fp;
	char buf[100];

	struct _Graph_node_A *PA;

	/***********************************************************
	  read in the trace file **********************************
	 ************************************************************/
	printf("Please input the trace file name:");
	scanf("%s", &tracefilename);

	fp = fopen(tracefilename, "r");
	if(fp==NULL){
		printf("Could not open the trace file!\n");
		exit(1);
	}
	printf("Has open the trace successfully!\n");

	while(fgets(buf,100,fp)){
		if(buf[0]=='%') continue;
		if(i==0){
			sscanf(buf, "%d%d%d", &max_id_A, &max_id_B, &max_weight);
			break;
		}
	}

	fclose(fp);

	cout << max_id_A << endl;
	PA = (struct _Graph_node_A *)malloc(sizeof(struct _Graph_node_A) * max_id_A);
	if(PA == NULL)
		cout << "Allocate memory for PA failed" <<endl;
	else
		cout << "Allocate memory for PA successfully" << endl;

	// initialize A and B
	for(i=0;i<max_id_A;i++)
		for(j=0;j<N;j++)
			PA[i].adj[j] = 0;

	cout << "Initialized PA successfully!" << endl;

	fp = fopen(tracefilename, "r");
	if(fp==NULL){
		printf("Could not open the trace file!\n");
		exit(1);
	}
	printf("Has open the trace sucessfully!\n");

	int index[max_id_A];
	for(i=0;i<max_id_A;i++)
		index[i]  = 0;
	cout << "Initilized index successfully" << endl;

	i=0;

	//read in the input file and build the graph
	while(fgets(buf,100,fp)){
		if(buf[0]=='%')
			continue;
		if(i==0){
			sscanf(buf,"%d%d%d", &max_id_A, &max_id_B, &max_weight);
			cout << max_id_A << " " << max_id_B << " " << max_weight << endl;
			i++;
		}else{
			sscanf(buf, "%d%d", &A, &B);
			PA[A-1].adj[index[A-1]++] = B;
		}
	}

	fclose(fp);

	// Naive Graph collaboratie Filtering
	int v_start;// the start vertive for search
	int k;// the k is the number of relational vertices needed to be found

	// read in the start vertex, and value k for ralation depth
	std::cout << "Please input the start vertice as an interger number less than " << max_id_A << " :";
	scanf("%d", &v_start);
	while(v_start<0 || v_start > max_id_A) {
		cout << "Please input a valid start vertice less than " << max_id_A << " :";
		scanf("%d", &v_start);
	}

	std::cout << "Please input the value of k as an interger number less than " << max_id_A << " :";
	scanf("%d", &k);
	while(k < 0 || k > max_id_A) {
		cout << "Please input a valid k less than " << max_id_A << " :";
		scanf("%d", &k);
	}


	int* relation_count = (int*)malloc((k+1)*max_id_A*sizeof(int));
	int* index_vertice = (int*)malloc((k+1)*max_id_A*sizeof(int));

	int* relation_count_h = (int*)malloc((k+1)*max_id_A*sizeof(int));
	int* index_vertice_h = (int*)malloc((k+1)*max_id_A*sizeof(int));

	//initialize the relation_count and vertice index
	for(i=0;i<k+1;i++)
		for(j=0;j<max_id_A;j++){
			relation_count[i * max_id_A + j] = 0;
			relation_count_h[i * max_id_A + j] = 0;
			index_vertice[i * max_id_A + j] = j + 1;
			index_vertice_h[i * max_id_A + j] = j + 1;
		}

	// executing the filtering algorithm on the host side
	clock_t st = clock();
	filterOnHost(PA, relation_count_h, index_vertice_h, v_start, k, max_id_A);
	st = clock() - st;
	printf("CPU execution time is %.5f\n", (float)st/CLOCKS_PER_SEC);

#if DEBUG
	cout << "Relation count on the host side is:" << endl;
	for(i=0;i<k+1;i++){
		for(j=0;j<max_id_A;j++){
			cout << relation_count_h[i*max_id_A+j] << " ";
		}
		cout << endl;
	}
#endif

	// executing the filtering algorithm on the device side
	st = clock();
	filterOnDevice(PA, relation_count, index_vertice, v_start, k, max_id_A);
	st = clock() - st;
	printf("GPU execution time is %.5f\n", (float)st/CLOCKS_PER_SEC);


	//sort the relation in a descending order
	for(i=1;i<k+1;i++) {
		selSort(&relation_count[i*max_id_A], &index_vertice[i*max_id_A], max_id_A);
		selSort(&relation_count_h[i*max_id_A], &index_vertice_h[i*max_id_A], max_id_A);
	}

#if DEBUG
	cout << "The " << k << " related vertices to each vertice are: " << endl;
	for(i=0;i<k+1;i++) {
		for(j=0;j<k;j++) {
			cout << index_vertice[i * max_id_A + j]	<< " ";
			cout << index_vertice_h[i * max_id_A + j] << " ";
		}
		cout << endl;
	}
#endif


	vector<int> final_result_d;// record the final result of collaborative filtering for visualization
	vector<int> final_result_h;
	vector<int>::iterator it;

	// copy the index_vertice to the final_result vector
	for(i=0;i<k+1;i++)
		for(j=0;j<k;j++) {
			final_result_d.push_back(index_vertice[i * max_id_A + j]);
			final_result_h.push_back(index_vertice_h[i * max_id_A + j]);
		}


	// sort the final_result vector in a desending order
	std::sort(final_result_d.begin(), final_result_d.end());
	std::sort(final_result_h.begin(), final_result_h.end());

	// remove the repeated vertices
	for(it=final_result_d.begin()+1, i=final_result_d.front();it!=final_result_d.end();) {
		if(i == *it)
			final_result_d.erase(it);
		else {
			i = *it;
			it++;
		}
	}

	for(it=final_result_h.begin()+1, i=final_result_h.front();it!=final_result_h.end();) {
		if(i == *it)
			final_result_h.erase(it);
		else {
			i = *it;
			it++;
		}
	}

	// compare the result from GPU with the result from CPU to test the correctness
	bool match = compare(final_result_h, final_result_d);
	if(match)
		cout << "Test passed ^^!" << endl;
	else
		cout << "Test failed !!" << endl;

	// output the final result
	cout << "The final Collaborative Filtering result is:" << endl;
	for(it=final_result_d.begin();it!=final_result_d.end();it++)
		cout << *it << " ";
	cout << endl;

	free(relation_count);
	free(relation_count_h);
	free(index_vertice);
	free(index_vertice_h);

	return 0;
}

void filterOnHost(struct _Graph_node_A *PA, int *relation_count, int *index_vertice, int v_start, int k, int max_id_A)
{
	int i, j, m, n;

	// naive collaborative filtering
	for(i=0;i<max_id_A;i++){
		if(i == v_start - 1)	continue;
		for(j=0;j<N && PA[i].adj[j] != 0;j++)
			for(m=0;m<N && PA[v_start-1].adj[m] != 0;m++){
				if(PA[v_start-1].adj[m] == PA[i].adj[j]) {
					relation_count[i]++;
					break;
				}
			}
	}

	// sort the relation in a descending order by selection sort algorithm
	selSort(relation_count, index_vertice, max_id_A);

	// full collaborative filtering
	for(i=0;i<k;i++)
		for(j=0;j<max_id_A;j++)
			for(m=0;m<N && PA[j].adj[m] != 0;m++)
				for(n=0;n<N && PA[index_vertice[i]-1].adj[n] != 0;n++) {
					if(PA[index_vertice[i]-1].adj[n] == PA[j].adj[m] && index_vertice[i] - 1 != j) {
						relation_count[(i+1)*max_id_A + j] ++;
						break;
					}
				}

}

void filterOnDevice(struct _Graph_node_A *PA, int *relation_count, int *index_vertice, int v_start, int k, int max_id_A) {
	struct _Graph_node_A *PA_d;
	int *relation_count_d;
	int *max_id_A_d;
	int *index_vertice_d;
	int *v_start_d;
	int *k_d;
	int size = (k + 1) * max_id_A * sizeof(int);

	cudaMalloc((void**)&(PA_d), sizeof(struct _Graph_node_A) * max_id_A);
	cudaMalloc((void**)&(relation_count_d), size);
	cudaMalloc((void**)&(max_id_A_d), sizeof(int));
	cudaMalloc((void**)&(v_start_d), sizeof(int));
	cudaMalloc((void**)&(k_d), sizeof(int));
	cudaMalloc((void**)&(index_vertice_d), max_id_A * sizeof(int));

	cudaMemcpy(PA_d, PA, sizeof(struct _Graph_node_A) * max_id_A, cudaMemcpyHostToDevice);
	cudaMemcpy(relation_count_d, relation_count, size, cudaMemcpyHostToDevice);
	cudaMemcpy(max_id_A_d, &max_id_A, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(v_start_d, &v_start, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(k_d, &k, sizeof(int), cudaMemcpyHostToDevice);

	// naive collaborative filtering
	dim3 dimGrid0(1, 1);
	dim3 dimBlock0(max_id_A, 1);
	naiveFilterKernel<<<dimGrid0, dimBlock0>>>(PA_d, relation_count_d, v_start_d, max_id_A_d);

	cudaMemcpy(relation_count, relation_count_d, max_id_A * sizeof(int), cudaMemcpyDeviceToHost);

#if DEBUG
	cout << "The relation count is:" << endl;
	for(int i=0; i<k+1;i++) {
		for(int j=0;j<max_id_A;j++){
			cout << relation_count[i*max_id_A+j] << " ";
		}
		cout << endl;
	}
#endif

	selSort(relation_count, index_vertice, max_id_A);

#if DEBUG
	cout << "The " << k << " related vertice to " << v_start << " is:" << endl;
	for(int i=0;i<k;i++)
		cout << index_vertice[i] << " ";
	cout << endl;
#endif

	cudaMemcpy(index_vertice_d, index_vertice, max_id_A * sizeof(int), cudaMemcpyHostToDevice);

	// full collaborative filtering
	dim3 dimGrid(1, k);
	dim3 dimBlock(max_id_A, 1);
	fullFilterKernel<<<dimGrid, dimBlock>>>(PA_d, relation_count_d, index_vertice_d, k_d, max_id_A_d);

	cudaMemcpy(relation_count, relation_count_d, size, cudaMemcpyDeviceToHost);

#if DEBUG
	cout << "The relation count is:" << endl;
	for(int i=0; i<k+1;i++) {
		for(int j=0;j<max_id_A;j++){
			cout << relation_count[i*max_id_A+j] << " ";
		}
		cout << endl;
	}
#endif

	cudaFree(PA_d);
	PA_d = NULL;
	cudaFree(relation_count_d);
	relation_count_d = NULL;
	cudaFree(max_id_A_d);
	max_id_A_d = NULL;
	cudaFree(k_d);
	k_d = NULL;
	cudaFree(index_vertice_d);
	index_vertice_d = NULL;
}

// Selection Sort
void selSort(int s[], int index[], int length)
{
	int i, j, maxPos;
	for(i=0;i<length-1;i++) {
		maxPos = i;
		for(j=i+1;j<length;j++)
			if(s[j] > s[maxPos])
				maxPos = j;
		if(i != maxPos) {
			swap(s[i], s[maxPos]);
			swap(index[i], index[maxPos]);
		}
	}
}

bool compare(vector<int> final_result_h, vector<int> final_result_d)
{
	vector<int>::iterator it_h;
	vector<int>::iterator it_d;

	for(it_h=final_result_h.begin(), it_d=final_result_d.begin(); it_h != final_result_h.end() && it_d != final_result_d.end();it_h++, it_d++) {
		if((*it_h) != (*it_d))
			return false;
	}

	return true;
}
