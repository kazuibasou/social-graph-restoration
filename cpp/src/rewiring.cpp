#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include "basic_function.h"
#include "graph.h"
#include "rewiring.h"

bool two_K_rewirable(const int &i, const int &j, const int &a, const int &b, const int &k_1, const int &k_2, const int &k_3, const int &k_4){
	if(i == j){
		return false;
	}
	else if(i == a){
		return false;
	}
	else if(i == b){
		return false;
	}
	else if(j == a){
		return false;
	}
	else if(j == b){
		return false;
	}
	else if(a == b){
		return false;
	}
	else if(k_1 == k_3){
		return true;
	}
	else if(k_1 == k_4){
		return true;
	}
	else if(k_2 == k_3){
		return true;
	}
	else if(k_2 == k_4){
		return true;
	}
	
	return false;
}

int calculate_num_tri_to_add(const Graph &G, const std::vector<int> &degrees, const int &i, const int &j, const int &b, std::vector<int> &num_tri_to_add){
	int t_minus, t_plus;

	for(int k: G.nlist[i]){
		if(degrees[k] <= 1 || i == k){continue;}

		if(j != k && degrees[j] > 1){
			t_minus = std::count(G.nlist[j].begin(), G.nlist[j].end(), k);
			num_tri_to_add[degrees[i]] -= t_minus;
			num_tri_to_add[degrees[j]] -= t_minus;
			num_tri_to_add[degrees[k]] -= t_minus;
		}

		if(b != k && degrees[b] > 1){
			t_plus = std::count(G.nlist[b].begin(), G.nlist[b].end(), k);
			num_tri_to_add[degrees[i]] += t_plus;
			num_tri_to_add[degrees[k]] += t_plus;
			num_tri_to_add[degrees[b]] += t_plus;
		}
	}

	return 0;
}

int rewiring_random_edge_with_preserving_two_K(const Graph &G, const std::vector<int> &degrees, const std::vector<std::pair<int, int>> &rewirable_edges, std::mt19937 &engine, std::uniform_int_distribution<int> &dist, std::vector<int> &num_tri_to_add, int &i_e1, int &i_e2, int &rewiring_case){
	i_e1 = dist(engine);
	i_e2 = dist(engine);
	std::pair<int, int> e1 = rewirable_edges[i_e1];
	std::pair<int, int> e2 = rewirable_edges[i_e2];

	int i = e1.first;
	int j = e1.second;
	int a = e2.first;
	int b = e2.second;

	int k_1 = degrees[i];
	int k_2 = degrees[j];
	int k_3 = degrees[a];
	int k_4 = degrees[b];

	while(! two_K_rewirable(i, j, a, b, k_1, k_2, k_3, k_4)){
		i_e1 = dist(engine);
		i_e2 = dist(engine);
		e1 = rewirable_edges[i_e1];
		e2 = rewirable_edges[i_e2];

		i = e1.first;
		j = e1.second;
		a = e2.first;
		b = e2.second;

		k_1 = degrees[i];
		k_2 = degrees[j];
		k_3 = degrees[a];
		k_4 = degrees[b];
	}

	int t_minus,c;
	if(k_1 == k_3 || k_2 == k_4){

		calculate_num_tri_to_add(G, degrees, i, j, b, num_tri_to_add);
		calculate_num_tri_to_add(G, degrees, a, b, j, num_tri_to_add);

		if(degrees[j] > 1 && degrees[b] > 1){
			t_minus = std::count(G.nlist[j].begin(), G.nlist[j].end(), b);
			num_tri_to_add[degrees[i]] -= t_minus;
			num_tri_to_add[degrees[j]] -= 2*t_minus;
			num_tri_to_add[degrees[a]] -= t_minus;
			num_tri_to_add[degrees[b]] -= 2*t_minus;
		}
		
		if(degrees[i] > 1 && degrees[a] > 1){
			c = std::count(G.nlist[a].begin(), G.nlist[a].end(), i);
			if(degrees[j] > 1){
				num_tri_to_add[degrees[i]] -= c;
				num_tri_to_add[degrees[a]] -= c;
				num_tri_to_add[degrees[j]] -= c;
			}
			if(degrees[b] > 1){
				num_tri_to_add[degrees[i]] -= c;
				num_tri_to_add[degrees[a]] -= c;
				num_tri_to_add[degrees[b]] -= c;
			}
		}

		rewiring_case = 0;
	}
	else if(k_1 == k_4 || k_2 == k_3){
		
		calculate_num_tri_to_add(G, degrees, i, j, a, num_tri_to_add);
		calculate_num_tri_to_add(G, degrees, b, a, j, num_tri_to_add);

		if(degrees[j] > 1 && degrees[a] > 1){
			t_minus = std::count(G.nlist[j].begin(), G.nlist[j].end(), a);
			num_tri_to_add[degrees[i]] -= t_minus;
			num_tri_to_add[degrees[j]] -= 2*t_minus;
			num_tri_to_add[degrees[a]] -= 2*t_minus;
			num_tri_to_add[degrees[b]] -= t_minus;
		}

		if(degrees[i] > 1 && degrees[b] > 1){

			c = std::count(G.nlist[b].begin(), G.nlist[b].end(), i);

			if(degrees[j] > 1){
				num_tri_to_add[degrees[i]] -= c;
				num_tri_to_add[degrees[b]] -= c;
				num_tri_to_add[degrees[j]] -= c;
			}
			
			if(degrees[a] > 1){
				num_tri_to_add[degrees[i]] -= c;
				num_tri_to_add[degrees[b]] -= c;
				num_tri_to_add[degrees[a]] -= c;
			}
		}

		rewiring_case = 1;
	}
	else{
		printf("Error: Selected edges destroy the present joint degree matrix.\n");
		exit(0);
	}

	return 0;
}

int calc_L1_distance(const std::vector<double> &target_ddcc, const std::vector<double> &current_ddcc, double &dist, double &norm){
    int target_s = int(target_ddcc.size());
    int current_s = int(current_ddcc.size());
    int s = std::max(target_s, current_s);
    dist = 0.0;
    norm = 0.0;

    for(int i=0; i<s; ++i){
    	if(i<target_s){
    		norm += target_ddcc[i];
    	}

    	if(i<target_s && i<current_s){
    		dist += std::fabs(target_ddcc[i]-current_ddcc[i]);
    	}
    	else if(i>=target_s && i<current_s){
    		dist += current_ddcc[i];
    	}
    	else if(i<target_s && i>=current_s){
    		dist += target_ddcc[i];
    	}
    	else{
    		printf("Error: Out of range.\n");
    		exit(0);
    	}
    }

    return 0;
}

int Rewiring::targeting_two_five_K(Graph &genG, const std::vector<double> &target_ddcc, std::vector<std::pair<int, int>> rewirable_edges){
	int i, d, id;
	int d_size = std::max(genG.maxd, int(target_ddcc.size())) + 1;

	std::vector<int> degrees(genG.N, 0);
	std::vector<int> V_d(d_size, 0);
	for(i=0; i<int(genG.N); ++i){
		d = int(genG.nlist[i].size());
		degrees[i] = d;
		V_d[d] += 1;
	}

	std::vector<double> coeff(d_size, 0.0);
	for(d=2; d<int(coeff.size()); ++d){
		if(V_d[d] == 0){continue;}
		coeff[d] = double(2)/(d*(d-1));
		coeff[d] = double(coeff[d])/V_d[d];
	}

	genG.calc_clustering();
	std::vector<double> current_ddcc(d_size, 0.0);
	for(id=0; id<int(genG.ddcc.entries.size()); ++id){
		d = genG.ddcc.id_to_keys[id];
		current_ddcc[d] = genG.ddcc.get(d);
	}

	double D, rewired_D, delta_D, norm;
	calc_L1_distance(target_ddcc, current_ddcc, D, norm);
	std::vector<double> rewired_ddcc(d_size, 0.0);
	std::vector<int> num_tri_to_add(d_size, 0);

	int i_e1, i_e2, rewiring_case, tmp;
	int R = R_C*int(rewirable_edges.size());
	int r;

	std::random_device seed_gen;
	std::mt19937 engine(seed_gen());
	std::uniform_int_distribution<int> dist(0, int(rewirable_edges.size())-1);

	//printf("Initial error: %lf\n", double(D)/norm);

	std::chrono::system_clock::time_point  start, end;
	start = std::chrono::system_clock::now();

	for(r=0; r<R; ++r){
		rewired_ddcc = std::vector<double>(current_ddcc);
		num_tri_to_add = std::vector<int>(d_size, 0);
		rewired_D = D;

		rewiring_random_edge_with_preserving_two_K(genG, degrees, rewirable_edges, engine, dist, num_tri_to_add, i_e1, i_e2, rewiring_case);

		for(d=2; d<int(num_tri_to_add.size()); ++d){
			if(num_tri_to_add[d] == 0){continue;}
			rewired_ddcc[d] += double(num_tri_to_add[d]*coeff[d]);
			rewired_D += std::fabs(target_ddcc[d] - rewired_ddcc[d]) - std::fabs(target_ddcc[d] - current_ddcc[d]);
		}

		delta_D = rewired_D - D;

		if(delta_D < 0){

			if(rewiring_case == 0){
				genG.remove_edge(rewirable_edges[i_e1].first, rewirable_edges[i_e1].second);
				genG.add_edge(rewirable_edges[i_e1].first, rewirable_edges[i_e2].second);
				genG.remove_edge(rewirable_edges[i_e2].first, rewirable_edges[i_e2].second);
				genG.add_edge(rewirable_edges[i_e1].second, rewirable_edges[i_e2].first);

				tmp = rewirable_edges[i_e1].second;
				rewirable_edges[i_e1].second = rewirable_edges[i_e2].second;
				rewirable_edges[i_e2].second = tmp;
			}
			else{
				genG.remove_edge(rewirable_edges[i_e1].first, rewirable_edges[i_e1].second);
				genG.add_edge(rewirable_edges[i_e1].first, rewirable_edges[i_e2].first);
				genG.remove_edge(rewirable_edges[i_e2].first, rewirable_edges[i_e2].second);
				genG.add_edge(rewirable_edges[i_e1].second, rewirable_edges[i_e2].second);

				tmp = rewirable_edges[i_e1].second;
				rewirable_edges[i_e1].second = rewirable_edges[i_e2].first;
				rewirable_edges[i_e2].first = tmp;
			}

			current_ddcc = std::vector<double>(rewired_ddcc);
			D = rewired_D;
		}
	}

	end = std::chrono::system_clock::now();
	//double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	//printf("Final error: %lf, Final acceptance rate: %lf\n", double(D)/norm, A);
	//printf("%lf [milli sec] %d [rewiring], %lf [millisec/rewiring]\n", elapsed, R*step, double(elapsed)/(R*step));

	return 0;
}