#ifndef REWIRING_H
#define REWIRING_H

class Rewiring{
	public:
		//parameter
	  	const int R_C = 500;

		Rewiring(){}
		~Rewiring(){}

		int targeting_two_five_K(Graph &genG, const std::vector<double> &target_ddcc, std::vector<std::pair<int, int>> rewirable_edges);
};

#endif