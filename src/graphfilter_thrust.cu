#include <iostream>
#include <utility>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>

//#include <cutil.h>

// thrust graph library
#include <thrust/graph/graph_traits.hpp>
#include <thrust/graph/adjacency_list.hpp>
#include <thrust/graph/random.hpp>
#include <thrust/graph/breadth_first_search.hpp>i
#include <thrust/graph/simple_io.hpp>
#include <thrust/graph/graph_selectors.hpp>
#include <thrust/graph/graph_concepts.hpp>
#include <thrust/graph/properties.hpp>
#include <thrust/graph/visitors.hpp>
#include "filter.hpp"

#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/reduce_by_keyvalue.hpp>

using namespace std;
using namespace thrust;

// function to sort the relations to the start vertice in a descending order
void selSort(int s[], int index[], int length);

static const char* DAT_FILE_NAME = "graph.dat";// file to store the graph struct

int main(int arcg, char** argv)
{
	typedef adjacency_list<directedS, disallow_parallelS> graph_t;
	//typedef adjacency_list<bidirectionalS, unparallelS, no_property, no_property, no_property> graph_t;
	typedef std::pair<int, int> Edge;
	//vector<Edge> edgeVec;
	std::size_t num_vertices;//total num of vertices in party one and party two
	char tracefilename[30];
	int A, B; //fot the nodes in party A and party B
	double Eweight;
	int max_id_A;
	int max_id_B;
	int max_weight;
	int i=0, j=0;
	int m=0, n=0;
	FILE *fp;
	char buf[100];

	//graph_t g(59);

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
		}
	}

	fclose(fp);

	num_vertices = max_id_A + max_id_B;
	graph_t g(num_vertices);

	fp = fopen(tracefilename, "r");
        if(fp==NULL){
                printf("Could not open the trace file!\n");
                exit(1);
        }
	printf("Has open the trace sucessfully!\n");
	while(fgets(buf,100,fp)){
                if(buf[0]=='%') continue;
                if(i==0){
                        sscanf(buf,"%d%d%d", &max_id_A, &max_id_B, &max_weight);
                        i++;
                }else{
                        sscanf(buf, "%d%d", &A, &B);
			thrust::add_edge(A, B+max_id_A, g);
			thrust::add_edge(B+max_id_A, A, g);
                        //edgeVec.push_back(Edge(A,B+max_id_A));
                }

        }

	fclose(fp);


	//declare and add the edges to the graph object
	//std::cout << g << endl;

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

	for(i=0;i<k+1;i++)
		for(j=0;j<max_id_A;j++){
			relation_count[i * max_id_A + j] = 0;
			index_vertice[i * max_id_A + j] = j + 1;
		}

	typedef graph_traits<graph_t> GraphTraits;
	typedef GraphTraits::vertex_descriptor Vertex;
	//typedef GraphTraits::vertex_iterator;
	typedef GraphTraits::edge_descriptor edge_descriptor;
	typedef GraphTraits::edge_iterator edge_iterator;


	GraphTraits::out_edge_iterator in_i, in_end; // in edges iterator for a vertice
	GraphTraits::adjacency_iterator ai, ai_end; // adjacent vertices iterator

	// get the property map for vertex indices
	typedef property_map<graph_t, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	//tie(ai, ai_end) = adjacent_vertices(v_start, g);
	//filter(g, relation_count, max_id_A, k, v_start);

	//thrust::host
	for(tie(ai, ai_end)= adjacent_vertices(v_start, g); ai != ai_end; ++ai)// get all the vertices adjacent to v_start
		for(tie(in_i, in_end) = out_edges(*ai, g); in_i != in_end; ++in_i) {
		edge_descriptor e = *in_i;
		Vertex src = target(e, g);
		if(src != v_start)
			relation_count[index[src]-1]++;
		}

		// sort the relation in a descending order by selection sort algorithm
		selSort(relation_count, index_vertice, max_id_A);

		std::cout << "the " << k << " relevant vertices to " << v_start << " is: " << std::endl;
		for(i=0;i<k;i++)
			std::cout << index_vertice[i] << " ";
		std::cout << std::endl;

	/***********************************************************************************
	  *****	Collaborative Filtering for visualization **********************************
	***********************************************************************************/
	for(i=0,j=1;i<k;i++,j++)
		for(tie(ai,ai_end) = adjacent_vertices(index_vertice[i], g); ai != ai_end; ++ai)
			for(tie(in_i, in_end) = out_edges(*ai, g); in_i != in_end; ++in_i) {
				edge_descriptor e = *in_i;
				Vertex src = target(e, g);
				if(src != index_vertice[i])
					relation_count[j*max_id_A + index[src] - 1]++;
			}


	// sort the relation in a descending order
	for(i=1;i<k+1;i++)
		selSort(&relation_count[i*max_id_A], &index_vertice[i*max_id_A], max_id_A);

	cout << "The " << k << " related vertices to each vertice are: " << endl;
	for(i=0;i<k+1;i++) {
		for(j=0;j<k;j++) {
			cout << index_vertice[i * max_id_A + j]	<< " ";
		}
		cout << endl;
	}


	vector<int> final_result;// record the final result of collaborative filtering for visualization
	vector<int>::iterator it;

	// copy the index_vertice to the final_result vector
	for(i=0;i<k+1;i++)
		for(j=0;j<k;j++)
			final_result.push_back(index_vertice[i * max_id_A + j]);

	// sort the final_result vector in a desending order
	std::sort(final_result.begin(), final_result.end());

	// remove the repeated vertices
	for(it=final_result.begin()+1, i=final_result.front();it!=final_result.end();) {
		if(i == *it)
			final_result.erase(it);
		else {
			i = *it;
			it++;
		}
	}

	// output the final result
	cout << "The final Collaborative Filtering result is:" << endl;
	for(it=final_result.begin();it!=final_result.end();it++)
		cout << *it << endl;

	free(relation_count);
	free(index_vertice);

	return 0;
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
