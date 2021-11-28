#ifndef SAMPLING_H
#define SAMPLING_H

class SampledData{
	public:
		int v;
		std::vector<int> nlist;
		SampledData(const int targetnode);
		~SampledData();
};

SampledData query(const Graph G, const int targetnode);

int select_seed(const Graph G);

int random_walk(const Graph G, const int samplesize,const int startnode, std::vector<SampledData> &samplinglist);

#endif