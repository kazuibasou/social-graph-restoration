#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string.h>
#include <random>
#include <vector>
#include <deque>
#include <stack>
#include <queue>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "basic_function.h"
#include "graph.h"

Graph::Graph(){
	N = 0;
	M = 0;
	maxd = 0;
	n_q = 0;
	n_v = 0;
	ddcc.clear();
}

Graph::~Graph(){

	std::vector<std::vector<int>>().swap(nlist);
}

int Graph::read_graph(const char *dir, const char *readgraph){

	int s = int(sizeof(readgraph))/sizeof(char);
	if(s > 100){
		printf("Error: Length of name of readfile must be less than 100.\n");
		exit(0);
	}

	FILE *rf;
	std::string rfpath = std::string(dir) + std::string(readgraph) + ".txt";
	rf = fopen(rfpath.c_str(), "r");
	if(rf == NULL) {
		printf("Error: Could not open file named %s.txt.\n",readgraph);
		exit(0);
	}

	int x,y;
	std::set<int> V;
	while(fscanf(rf, "%d %d", &x , &y) != EOF){
		V.insert(x);
		V.insert(y);
	}
	fclose(rf);

	N = int(V.size());
	nlist = std::vector<std::vector<int>>(N,std::vector<int>());
	rf = fopen(rfpath.c_str(), "r");
	while(fscanf(rf, "%d %d", &x , &y) != EOF){
		if(x < 0 || N <= x){
			printf("Error: Node index must be between 0 and N-1.\n");
			exit(0);
		}
		if(y < 0 || N <= y){
			printf("Error: Node indesx must be between 0 and N-1.\n");
			exit(0);
		}
		nlist[x].push_back(y);
		nlist[y].push_back(x);
	}
	fclose(rf);

	M = 0;
	maxd = 0;
	for(int i=0;i<N;++i){
		M += int(nlist[i].size());
		std::sort(nlist[i].begin(),nlist[i].end());
		if(int(nlist[i].size()) > maxd){
			maxd = nlist[i].size();
		}
	}
	M = int(M)/2;

	printf("Network data: %s\n", readgraph);
	printf("Number of nodes: %d, Number of edges: %d\n", N, M);

	return 0;
}

int Graph::add_edge(const int v,const int w){
	//Add an edge between v and w.

	/*
	if(int(std::max(v,w)+1) > int(nlist.size())){
		nlist.resize(std::max(v,w)+1, std::vector<int>());
	}
	*/

	nlist[v].push_back(w);
	nlist[w].push_back(v);
	
	return 0;
}

int Graph::remove_edge(const int v,const int w){
	//Remove an edge between v and w.
	
	auto itr1 = std::find(nlist[v].begin(), nlist[v].end(), w);

	/*
	if(itr1 == nlist[v].end()){
		printf("Error: node w is not included in neighbors of node v.\n");
		exit(0);
	}
	*/
	
	nlist[v].erase(itr1);

	
	auto itr2 = std::find(nlist[w].begin(), nlist[w].end(), v);
	
	/*
	if(itr2 == nlist[w].end()){
		printf("Error: node v is not included in neighbors of node w.\n");
		exit(0);
	}
	*/
	
	nlist[w].erase(itr2);
	
	return 0;
}

int Graph::calc_clustering(){
	int i, j, k, d;
	Vector<int> V_d;
	Vector<double> sum_cc_d;

	for(i=0; i<N; ++i){
		double cc = 0.0;
		d = nlist[i].size();
		V_d.add(d, 1);

		if(d == 0 || d == 1){continue;}

		for(int ij=0; ij<d-1; ++ij){
			j = nlist[i][ij];
			for(int ik=ij+1; ik<d; ++ik){
				k = nlist[i][ik];
				if(i != j && j != k && k != i){
					cc += 2*std::count(nlist[j].begin(), nlist[j].end(), k);
				}
			}
		}

		cc = double(cc)/(d*(d-1));
		sum_cc_d.add(d, cc);
	}

	ddcc.clear();
	for(i=0; i<int(V_d.entries.size()); ++i){
		d = V_d.id_to_keys[i];
		if(V_d.get(d) > 0){
			ddcc.set(d, double(sum_cc_d.get(d))/(V_d.get(d)));
		}
	}
	
	return 0;
}
