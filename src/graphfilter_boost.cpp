#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <vector>
// boost library
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/bipartite.hpp>

// arhivers
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace boost;
using namespace std;
//#define DEBUG 0;

void selSort(int s[], int index[], int length);// function to sort the relations to start vertice in a descending order


static const char* DAT_FILE_NAME = "graph.dat";// file to store the graph struct
  
int main(int,char*[])
{
 
	// create a typedef for the Graph type
 	typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
	typedef std::pair<int, int> Edge;
	vector<Edge> edgeVec;
    	int num_vertices;// total number of vertices in party one and party two
	char tracefilename[30];
	int A, B;//for the nodes in party A and party B
	double Eweight;
	int max_id_A;
	int max_id_B;
	int max_weight;
	int i=0, j=0;
	int m=0, n=0;
	FILE *fp;
	char buf[100];	

	// *********************************************************************************
	// read in the trace file **********************************************************
	printf("Please input trace file name:");
  	scanf("%s",&tracefilename);	
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
			edgeVec.push_back(Edge(A,B+max_id_A));
		}
 
        }
	fclose(fp);
	
	num_vertices = max_id_A + max_id_B;// the total number of vertices in both the two parties
	
 	// declare and add the edges to the graph object
 	Graph g(edgeVec.begin(), edgeVec.end(), num_vertices);
	

	// test if the Graph is a bipartite Graph
	if(!is_bipartite(g)){
		std::cerr << "The graph is not a bipartite graph !!" << std::endl;
		return EXIT_FAILURE;
	}	


	// graph serialization
	// serialize and save graph
	std::ofstream ofs(DAT_FILE_NAME, std::ios::out | std::ios::binary);
	if(!ofs.is_open()) {
		std::cerr << "Can't open " << DAT_FILE_NAME << " file." << std::endl;
		return EXIT_FAILURE;
	}
	boost::archive::binary_oarchive oa(ofs);
	oa << g;
	ofs.close();


	//reload the grap from the archive file
	std::ifstream ifs(DAT_FILE_NAME, std::ios::in | std::ios::binary);
        if(!ifs.is_open()) {
                std::cerr << "Can't open " << DAT_FILE_NAME << " file." << std::endl;
                return EXIT_FAILURE;
        }
        boost::archive::binary_iarchive ia(ifs);
	Graph g1;
	ia >> g1;


	// Naive Graph collaborative Filtering
	int v_start;// the start vertice for search
	int k;// the k is the number of relational vertices need to be found
	vector<int> final_result;// record the final result of collaborative filtering for visualization
	vector<int>::iterator it;
	
	std::cout << "Please input the start vertice as an interger number less than " << max_id_A << " :";
	scanf("%d", &v_start);
	while(v_start<0 || v_start > max_id_A){
		cout << "Please input a valid start vertice less than " << max_id_A << " :";
		scanf("%d", &v_start);
	}
	
	std::cout << "Please input the value of k as an interger number less than " << max_id_A << " :";
	scanf("%d", &k);
	while(k < 0 || k > max_id_A){
                cout << "Please input a valid K less than " << max_id_A << " :";
                scanf("%d", &k);
        }

	/* Define and initialize a two dimension dynamic matrix */
	int **relation_count;// To record the relation between the vertices
	int **index_vertice;
	relation_count = (int **)malloc(sizeof(int *) * (k+1));
	index_vertice = (int **)malloc(sizeof(int *) * (k+1));

	for(i=0;i<k+1;i++){
		relation_count[i] = (int *)malloc(sizeof(int) * max_id_A);
		index_vertice[i] = (int *)malloc(sizeof(int) * max_id_A);
		for(j=0;j<max_id_A;j++){
			relation_count[i][j] = 0;
			index_vertice[i][j] = j+1;
		}
			
	}

	
		
 
	typedef graph_traits<Graph> GraphTraits;
      	
	GraphTraits::in_edge_iterator in_i, in_end; // in edges iterator for a vertice
        GraphTraits::adjacency_iterator ai, ai_end;// adjacent vertices iterator
	GraphTraits::edge_descriptor e;	// edge
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	
	// get the property map for vertex indices
	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	clock_t st = clock();
	for(tie(ai,ai_end) = adjacent_vertices(v_start,g1); ai != ai_end; ++ai)	// get all the vertices adjacent to the start vertice
		for(tie(in_i,in_end) = in_edges(*ai,g1); in_i != in_end; ++in_i){ // if a vertice share a same nabor with the start vertice then the related counter increased by one
        		e = *in_i;
        		Vertex src = source(e, g1);
			if(src != v_start)		
				relation_count[0][index[src]-1]++;
		}
	st = clock() - st;
	printf("Latency of relation count for naive filtering is %.5f\n", (float)st/CLOCKS_PER_SEC);
 	
	st = clock();
	selSort(relation_count[0], index_vertice[0], max_id_A);// sort the relation in a descending order by selection sort algorithm
	st = clock() - st;
	printf("Latency of sorting for naive filtering is %.5f\n", (float)st/CLOCKS_PER_SEC);

#if 0	
	std::cout << "the " << k << " relevant vertices to " << v_start << " is: " << std::endl;
	for(i=0;i<k;i++)
		std::cout << index_vertice[0][i] << " ";
	std::cout << std::endl; 	
#endif	

	// ************************************************************************************************************************************
	// Collaborative Filtering for Visualization
	// **************************************************************************************************************************************
	/* calculate the relation betweent the start vertice and the other vertices */
	st = clock();
	for(i=0,j=1;i<k;i++,j++)
		for(tie(ai,ai_end) = adjacent_vertices(index_vertice[0][i],g1);ai != ai_end; ++ai)
			for(tie(in_i,in_end) = in_edges(*ai,g1); in_i != in_end; ++in_i){
				e = *in_i;
				Vertex src = source(e, g1);
				if(src != index_vertice[0][i])
					relation_count[j][index[src]-1]++;
			}
	st = clock() - st;
        printf("Latency of relation count for full filtering is %.5f\n", (float)st/CLOCKS_PER_SEC);

	/* sort the relation in a descending order */
	st = clock();		
	for(i=1;i<k+1;i++)
		selSort(relation_count[i], index_vertice[i], max_id_A);
	st = clock() - st;
	printf("Latency of sorting for full filtering is %.5f\n", (float)st/CLOCKS_PER_SEC);
	
#if 0
	cout << "The k related vertices to each vertice are:" << endl;
	for(i=0;i<k+1;i++){
		for(j=0;j<k;j++){
			cout << index_vertice[i][j] << " ";
		}
		cout << endl;
	}
#endif	
	/* copy the index_vertice to the final_result vector */
	for(i=0;i<k+1;i++)
	for(j=0;j<k;j++)
		final_result.push_back(index_vertice[i][j]);

	/* sort the final_result vector in a secending order  */
	st = clock();
	sort(final_result.begin(), final_result.end());

	/*  remove the repeated vertices  */
	for(it=final_result.begin()+1, i=final_result.front();it!=final_result.end();){
		if(i==*it)
			final_result.erase(it);
		else {
			i = *it;
			it++;
		}
	}		
	st = clock() - st;
	printf("Latency of repeated vertices removal for full filtering is %.5f\n", (float)st/CLOCKS_PER_SEC);	
#if 0	
	/* output the final result */
	cout << "The final Collaborative Filtering result is:" << endl;
	for(it=final_result.begin();it!=final_result.end();it++)
		cout << *it << endl;
#endif
	for(i=0;i<k+1;i++){
		free((void *)relation_count[i]);
		free((void *)index_vertice[i]);
	}
	
	return 0;
}



//Selection Sort
void selSort(int s[], int index[], int length)
{
	int i, j, maxPos;
	for(i=0;i<length-1;i++){
		maxPos = i;
		for(j=i+1;j<length;j++)
			if(s[j] > s[maxPos])
				maxPos = j;
		if(i != maxPos){
			swap(s[i],s[maxPos]);
			swap(index[i],index[maxPos]);
		}
	}
}
