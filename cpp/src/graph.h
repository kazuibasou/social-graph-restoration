#ifndef GRAPH_H
#define GRAPH_H

class Graph{
	public:
		int N; //number of nodes
		int M; //number of edges
		std::vector<std::vector<int>> nlist; //neighbor lists
		int maxd; //maximum degree
		int n_q; //number of queried nodes in a subgraph obtained through sampling
		int n_v; //number of visible nodes in a subgraph obtained through sampling
		Vector<double> ddcc; //degree-dependent clustering coefficient

		Graph();
		~Graph();

		int read_graph(const char *dir, const char *readgraph);
		int add_edge(const int v, const int w);
		int remove_edge(const int v, const int w);

		int calc_clustering();
};

#endif