#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <numeric>
#include "basic_function.h"
#include "graph.h"
#include "sampling.h"
#include "estimation.h"
#include "generation.h"
#include "rewiring.h"

int write_edge_list(const Graph &genG, const char *readgraph){

	const char *dir = "../../gen_graph/";
	const std::string file_path = std::string(dir) + std::string(readgraph) + ".txt";
	FILE *f = fopen(file_path.c_str(),"w");
	if(f == NULL) {
		printf("Could not open file named %s.\n", file_path.c_str());
		exit(0);
	}

	for(int i=0; i<genG.N; ++i){
		for(int j:genG.nlist[i]){
			if(j >= i){
				fprintf(f, "%d %d\n", i, j);
			}
		}
	}

	fclose(f);

	return 0;
}

int main(int argc,char *argv[]){
	if(argc != 3){
		printf("Please input following: ./main (name of dataset) (length of a simple random walk)\n");
		exit(0);
	}
	const char *dir_dataset = "../../data/";
	const char *readgraph = argv[1];
	const int samplesize = std::atoi(argv[2]);

	Graph G;
	G.read_graph(dir_dataset,readgraph);

	if(samplesize <= 0 || samplesize > G.N){
		printf("Error: The number of nodes to be queried must be not less than 0 and not more than %d.\n", G.N);
		exit(0);
	}

	srand((unsigned) time(NULL));

	int seed = select_seed(G);

	printf("The number of queried nodes: %d (The percentage of queried nodes: %lf)\n", samplesize, double(samplesize)/G.N);
	std::vector<SampledData> samplinglist;
	random_walk(G, samplesize, seed, samplinglist);
	
	Graph genG = graph_restoration_method(samplinglist);
	write_edge_list(genG, readgraph);

	return 0;
}