#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <numeric>
#include <random>
#include <algorithm>
#include "basic_function.h"
#include "graph.h"
#include "sampling.h"

SampledData::SampledData(const int targetnode){
	v = targetnode;
	nlist = std::vector<int>();
}

SampledData::~SampledData(){
	
	std::vector<int>().swap(nlist);
}

SampledData query(const Graph G, const int targetnode){
	SampledData data(targetnode);
	data.nlist = G.nlist[targetnode];

	return data;
}

int select_seed(const Graph G){

	return generate_rand(G.N);
}

//Simple random walk
int random_walk(const Graph G, const int samplesize, const int startnode, std::vector<SampledData> &samplinglist){
	
	if(samplesize <= 0 || samplesize > G.N){
		printf("Error: The number of nodes to be queried must be not less than 0 and not more than %d.\n", G.N);
		exit(0);
	}
	
	int v = startnode;
	int nextv,index;
	samplinglist.clear();
	std::vector<int> queried(G.N, 0);
	std::vector<int> queried_nodes;

	while(int(queried_nodes.size()) < samplesize){
		SampledData data = query(G, v);
		samplinglist.push_back(data);
		if(queried[v] == 0){
			queried[v] = 1;
			queried_nodes.push_back(v);
		}
		index = generate_rand(data.nlist.size());
		nextv = data.nlist[index];
		v = nextv;
	}

	if(int(queried_nodes.size()) != samplesize){
		printf("Error: RW did not successfully collect target number of samples.\n");
		exit(0);
	}

	return 0;
}