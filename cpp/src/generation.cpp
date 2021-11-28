#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <unordered_map>
#include <climits>
#include <algorithm>
#include <deque>
#include <cmath>
#include <limits>
#include <numeric>
#include "basic_function.h"
#include "graph.h"
#include "sampling.h"
#include "estimation.h"
#include "generation.h"
#include "rewiring.h"

int initialize_num_deg(const double est_n, const Vector<double> est_dd, Vector<int> &num_deg){
	//initialize numbers of nodes with each degree.

	num_deg = Vector<int>();
	int i, d, x;

	for(i=0; i<int(est_dd.entries.size()); ++i){
		d = est_dd.id_to_keys[i];
		x = my_round(est_dd.get(d)*est_n);
		num_deg.set(d, std::max(x, 1));
	}

	return 0;
}

int adjust_num_deg(const double est_n, const Vector<double> est_dd, Vector<int> &target_num_deg){
	//make sum of degrees an even number.

	int i, d, id;
	int sum_d = 0;
	double x, y, delta_e;

	for(i=0; i<int(target_num_deg.entries.size()); ++i){
		d = target_num_deg.id_to_keys[i];
		sum_d += d*target_num_deg.get(d);
	}

	if(sum_d % 2 == 0){return 0;}

	std::vector<std::pair<int, double>> degree_candidate;
	for(id=0; id<int(target_num_deg.entries.size()); ++id){
		d = target_num_deg.id_to_keys[id];

		if(d % 2 == 0){continue;}

		x = est_dd.get(d)*est_n;
		y = target_num_deg.get(d);

		if(x != 0){
			delta_e = double(std::fabs(x-y-1))/x - double(std::fabs(x-y))/x;
		}
		else{
			delta_e = INFINITY;
		}

		degree_candidate.push_back({d, delta_e});
	}

	if(int(degree_candidate.size()) > 0){
		d = min_key_with_smallest_value(degree_candidate);
		target_num_deg.add(d, 1);
	}
	else{
		target_num_deg.add(1, 1);
	}

	return 0;
}

bool compare_second(const std::pair<int, int> a, const std::pair<int, int> b){
	if(a.second != b.second){
		return a.second > b.second;
	}
	else{
		return a.first < b.first;
	}
}

int modify_num_deg(const Graph genG, const double est_n, const Vector<double> est_dd, 
	std::vector<int> &assigned_degree, Vector<int> &current_num_deg, Vector<int> &target_num_deg){

	int v, subG_d, id, k;
	current_num_deg = Vector<int>();
	for(v=0; v<genG.n_q; ++v){
		subG_d = int(genG.nlist[v].size());
		assigned_degree[v] = subG_d;
		current_num_deg.add(assigned_degree[v], 1);
	}

	for(id=0; id<int(current_num_deg.id_to_keys.size()); ++id){
		k = current_num_deg.id_to_keys[id];
		if(target_num_deg.get(k) < current_num_deg.get(k)){
			target_num_deg.set(k, current_num_deg.get(k));
		}
	}

	int subG_N = genG.n_q + genG.n_v;
	std::vector<std::pair<int,int>> visible_node_pairs;
	for(v=genG.n_q; v<subG_N; ++v){
		visible_node_pairs.push_back({v, int(genG.nlist[v].size())});
	}
	std::sort(visible_node_pairs.begin(), visible_node_pairs.end(), compare_second);

	for(auto visible_node_pair: visible_node_pairs){
		v = visible_node_pair.first;
		subG_d = visible_node_pair.second;

		std::vector<int> degree_candidates;
		for(id=0; id<int(target_num_deg.id_to_keys.size()); ++id){
			k = target_num_deg.id_to_keys[id];
			if(k >= subG_d && target_num_deg.get(k) > current_num_deg.get(k)){
				for(int i=0; i<int(target_num_deg.get(k) - current_num_deg.get(k)); ++i){
					degree_candidates.push_back(k);
				}
			}
		}

		if(int(degree_candidates.size()) > 0){
			assigned_degree[v] = random_choice(degree_candidates);
		}
		else{
			std::vector<std::pair<int, double>> degree_to_add_candidates;
			double x, y, delta_e;
			for(id=0; id<int(est_dd.entries.size()); ++id){
				k = est_dd.id_to_keys[id];
				if(k < subG_d){continue;}

				x = est_dd.get(k)*est_n;
				y = double(target_num_deg.get(k));
				if(x != 0){
					delta_e = double(std::fabs(x-y-1))/x - double(std::fabs(x-y))/x;
				}
				else{
					delta_e = INFINITY;
				}

				degree_to_add_candidates.push_back({k, delta_e});
			}

			if(int(degree_to_add_candidates.size()) > 0){
				assigned_degree[v] = min_key_with_smallest_value(degree_to_add_candidates);
			}
			else{
				assigned_degree[v] = subG_d;
			}
			target_num_deg.add(assigned_degree[v], 1);
		}

		current_num_deg.add(assigned_degree[v], 1);
	}

	adjust_num_deg(est_n, est_dd, target_num_deg);

	return 0;
}

int initialize_num_jnt_deg(const double est_n, const double est_aved, const Matrix<double> est_jdd, 
	Matrix<int> &target_num_jnt_deg){
	
	//initialize numbers of edges with joint degrees

	target_num_jnt_deg = Matrix<int>();
	int i, j, k, l, x;

	for(i=0; i<int(est_jdd.entries.size()); ++i){
		k = est_jdd.id_to_keys[i];
		for(j=0; j<int(est_jdd.entries[i].size()); ++j){
			l = est_jdd.id_to_keys[j];

			if(est_jdd.get(k, l) <= 0){continue;}

			x = my_round(est_n*est_aved*est_jdd.get(k, l));

			if(k != l){
				target_num_jnt_deg.set(k, l, std::max(x, 1));
			}
			else{
				if(x % 2 == 0){
					target_num_jnt_deg.set(k, l, std::max(x, 2));
				}
				else{
					if(std::fabs(est_n*est_aved*est_jdd.get(k, l) - x + 1) <= std::fabs(est_n*est_aved*est_jdd.get(k, l) - x - 1)){
						target_num_jnt_deg.set(k, l, std::max(x - 1, 2));
					}
					else{
						target_num_jnt_deg.set(k, l, std::max(x + 1, 2));
					}
				}
			}
		}
	}

	return 0;
}

int adjust_num_jnt_deg(const double est_n, const double est_aved, const Matrix<double> est_jdd, 
	const Matrix<int> &current_num_jnt_deg, Vector<int> &target_num_deg, Matrix<int> &target_num_jnt_deg){

	std::vector<int> degrees;
	target_num_deg.get_keys(degrees);
	if(std::find(degrees.begin(), degrees.end(), 1) == degrees.end()){
		degrees.push_back(1);
	}
	std::sort(degrees.begin(), degrees.end(), std::greater<int>());

	//make sum of m(k, l) for all l equal to k*n(k)
	int l, target_sum, current_sum, change;
	double x, y, delta_e;
	for(int k:degrees){
		target_sum = k*target_num_deg.get(k);
		current_sum = target_num_jnt_deg.sum_line(k);
		change = target_sum - current_sum;

		if(change == 0){continue;}

		std::vector<int> ls;
		for(int l:degrees){
			if(l <= k){
				ls.push_back(l);
			}
		}

		if(k == 1 && std::abs(target_sum - current_sum) % 2 != 0){
			target_num_deg.add(1, 1);
			target_sum += 1;
		}

		while(target_sum - current_sum != 0){

			if(target_sum > current_sum){
				std::vector<std::pair<int, double>> degree_candidate;
				for(int l:ls){

					if(current_sum == target_sum - 1 && l == k){continue;}

					x = est_jdd.get(k, l)*est_n*est_aved;
					y = double(target_num_jnt_deg.get(k, l));

					if(x == 0){
						delta_e = INFINITY;
					}
					else{
						if(l != k){
							delta_e = double(std::fabs(x - y - 1))/x - double(std::fabs(x - y))/x;
						}
						else{
							delta_e = double(std::fabs(x - y - 2))/x - double(std::fabs(x - y))/x;
						}
					}

					degree_candidate.push_back({l, delta_e});
				}

				l = random_key_with_smallest_value(degree_candidate);

				target_num_jnt_deg.add(k, l, 1);
				target_num_jnt_deg.add(l, k, 1);

				if(k != l){
					current_sum += 1;
				}
				else{
					current_sum += 2;
				}
			}
			else{
				std::vector<std::pair<int, double>> degree_candidate;
				for(int l:ls){
					if(target_num_jnt_deg.get(k, l) - current_num_jnt_deg.get(k, l) <= 0){continue;}

					if(current_sum == target_sum + 1 && l == k){continue;}

					x = est_jdd.get(k, l)*est_n*est_aved;
					y = double(target_num_jnt_deg.get(k, l));

					if(x == 0.0){
						delta_e = INFINITY;
					}
					else{
						if(l != k){
							delta_e = double(std::fabs(x - y + 1))/x - double(std::fabs(x - y))/x;
						}
						else{
							delta_e = double(std::fabs(x - y + 2))/x - double(std::fabs(x - y))/x;
						}
					}

					degree_candidate.push_back({l, delta_e});
				}

				if(int(degree_candidate.size()) > 0){
					l = random_key_with_smallest_value(degree_candidate);
					target_num_jnt_deg.subtract(k, l, 1);
					target_num_jnt_deg.subtract(l, k, 1);
					
					if(k != l){
						current_sum -= 1;
					}
					else{
						current_sum -= 2;
					}
				}
				else if(k > 1){
					target_sum += k;
					target_num_deg.add(k, 1);
				}
				else{
					target_sum += 2;
					target_num_deg.add(1, 2);
				}
			}
		}
	}

	return 0;
}

int modify_num_jnt_deg(const Graph genG, const double est_n, const double est_aved, const Matrix<double> est_jdd,
	const std::vector<int> assigned_degree, Vector<int> &target_num_deg, Matrix<int> &current_num_jnt_deg, 
	Matrix<int> &target_num_jnt_deg){

	std::vector<int> degrees;
	target_num_deg.get_keys(degrees);
	if(std::find(degrees.begin(), degrees.end(), 1) == degrees.end()){
		degrees.push_back(1);
	}

	int v, k1, k2, subG_N, i, j;
	current_num_jnt_deg = Matrix<int>();
	subG_N = genG.n_q + genG.n_v;
	for(v=0; v<subG_N; ++v){
		k1 = assigned_degree[v];
		for(int w:genG.nlist[v]){
			k2 = assigned_degree[w];
			current_num_jnt_deg.add(k1, k2, 1);
		}
	}

	int k3, k4;
	double x, y, delta_e;
	for(i=0; i<int(current_num_jnt_deg.entries.size()); ++i){
  		k1 = current_num_jnt_deg.id_to_keys[i];
  		for(j=0; j<int(current_num_jnt_deg.entries[i].size()); ++j){
  			k2 = current_num_jnt_deg.id_to_keys[j];

  			while(current_num_jnt_deg.get(k1, k2) > target_num_jnt_deg.get(k1, k2)){
  				target_num_jnt_deg.add(k1, k2, 1);
				target_num_jnt_deg.add(k2, k1, 1);

				std::vector<std::pair<int, double>> k3_candidates;
				for(int k3:degrees){
					if(k3 == k1){continue;}

					x = est_jdd.get(k2, k3)*est_n*est_aved;
					y = double(target_num_jnt_deg.get(k2, k3));

					if(y <= current_num_jnt_deg.get(k2, k3)){continue;}

					if(x == 0){
						delta_e = INFINITY;
					}
					else{
						if(k2 != k3){
							delta_e = double(std::fabs(x - y + 1))/x - double(std::fabs(x - y))/x;
						}
						else{
							delta_e = double(std::fabs(x - y + 2))/x - double(std::fabs(x - y))/x;
						}
					}

					k3_candidates.push_back({k3, delta_e});
				}

				k3 = 0;
				if(int(k3_candidates.size()) > 0){
					k3 = random_key_with_smallest_value(k3_candidates);
					target_num_jnt_deg.subtract(k2, k3, 1);
					target_num_jnt_deg.subtract(k3, k2, 1);
				}

				std::vector<std::pair<int, double>> k4_candidates;
				for(int k4:degrees){
					if(k4 == k2){continue;}

					x = est_jdd.get(k1, k4)*est_n*est_aved;
					y = double(target_num_jnt_deg.get(k1, k4));

					if(y <= current_num_jnt_deg.get(k1, k4)){continue;}

					if(x == 0){
						delta_e = INFINITY;
					}
					else{
						if(k1 != k4){
							delta_e = double(std::fabs(x - y + 1))/x - double(std::fabs(x - y))/x;
						}
						else{
							delta_e = double(std::fabs(x - y + 2))/x - double(std::fabs(x - y))/x;
						}
					}

					k4_candidates.push_back({k4, delta_e});
				}

				if(int(k4_candidates.size()) > 0){
					k4 = random_key_with_smallest_value(k4_candidates);
					target_num_jnt_deg.subtract(k4, k1, 1);
					target_num_jnt_deg.subtract(k1, k4, 1);

					if(k3 > 0){
						target_num_jnt_deg.add(k3, k4, 1);
						target_num_jnt_deg.add(k4, k3, 1);
					}
				}
  			}
  		}
  	}

  	adjust_num_jnt_deg(est_n, est_aved, est_jdd, current_num_jnt_deg, target_num_deg, target_num_jnt_deg);

  	return 0;
}

Graph graph_restoration_method(const std::vector<SampledData> samplinglist){

	Graph genG;

	//(1) construct subgraph
  	std::unordered_map<int,int> index;
  	std::vector<std::pair<int, int>> edges_to_add;
  	
	//(1-1) index queried nodes
	int i = 0;
	for(SampledData data:samplinglist){
		int v = data.v;
		if(index.find(v) == index.end()){
			index[v] = i;
			++i;
		}
	}
	genG.n_q = i;

	//(1-2) index visible nodes
	std::vector<int> sampled;
	for(SampledData data:samplinglist){
		int v = data.v;
		
		if(std::find(sampled.begin(), sampled.end(), v) != sampled.end()){continue;}
		
		sampled.push_back(v);

		for(int w:data.nlist){
			if(index.find(w) == index.end()){
				index[w] = i;
				++i;
			}
			if(std::find(sampled.begin(), sampled.end(), w) == sampled.end()){
				edges_to_add.push_back({index[v], index[w]});
			}
		}
	}
	genG.n_v = i - genG.n_q;

	//(1-3) construct subgraph induced from sampling list
	genG.nlist.resize(genG.n_q + genG.n_v);
	for(std::pair<int, int> e:edges_to_add){
		genG.add_edge(e.first, e.second);
	}

	//(1-4) If there are no visible nodes, return the subgraph.
	if(genG.n_v == 0){
		genG.N = genG.n_q + genG.n_v;
		genG.M = 0;
	  	genG.maxd = 0;
	  	for(int i=0;i<genG.N;++i){
	  		genG.M += genG.nlist[i].size();
	  		if(int(genG.nlist[i].size()) > genG.maxd){
	  			genG.maxd = genG.nlist[i].size();
	  		}
	  	}
	  	genG.M = int(genG.M)/2;

	  	return genG;
	}

	//(2) make target numbers of nodes with degree k
	double est_n;
	size_estimator(samplinglist, est_n);
	Vector<double> est_dd;
	degree_distribution_estimator(samplinglist, est_dd);

	Vector<int> target_num_deg;
	initialize_num_deg(est_n, est_dd, target_num_deg);
	adjust_num_deg(est_n, est_dd, target_num_deg);

	std::vector<int> assigned_degree(genG.n_q + genG.n_v, 0);
	Vector<int> current_num_deg;
	modify_num_deg(genG, est_n, est_dd, assigned_degree, current_num_deg, target_num_deg);

	//(3) make target numbers of edges between nodes with degree k and nodes with degree k'
	double est_aved;
	average_degree_estimator(samplinglist, est_aved);
	Matrix<double> est_jdd;
	jdd_estimator_hybrid(samplinglist, est_n, est_aved, est_jdd);

	Matrix<int> target_num_jnt_deg, current_num_jnt_deg;
	initialize_num_jnt_deg(est_n, est_aved, est_jdd, target_num_jnt_deg);
	adjust_num_jnt_deg(est_n, est_aved, est_jdd, current_num_jnt_deg, target_num_deg, target_num_jnt_deg);
	modify_num_jnt_deg(genG, est_n, est_aved, est_jdd, assigned_degree, target_num_deg, current_num_jnt_deg, target_num_jnt_deg);

	//(4) stub-based construction algorithm
	//(4-1) decide target number of nodes in a generated graph
	int target_n = 0;
	for(int id=0; id<int(target_num_deg.entries.size()); ++id){
		int d = target_num_deg.id_to_keys[id];
		target_n += target_num_deg.get(d);
	}

	genG.N = target_n;
	genG.nlist.resize(target_n);
	assigned_degree.resize(target_n);

	//(4-2) assign degree to added node
	std::vector<int> deg_seq;
	for(int id=0; id<int(target_num_deg.entries.size()); ++id){
		int d = target_num_deg.id_to_keys[id];
		for(int n_d=0; n_d<target_num_deg.get(d)-current_num_deg.get(d); ++n_d){
			deg_seq.push_back(d);
		}
	}

	std::random_device seed_gen;
  	std::mt19937 engine(seed_gen());
  	std::shuffle(deg_seq.begin(), deg_seq.end(), engine);

  	for(int i=0; i<int(deg_seq.size()); ++i){
  		int v = i+genG.n_q+genG.n_v;
  		assigned_degree[v] = deg_seq[i];
  		current_num_deg.add(assigned_degree[v], 1);
  	}

  	//(4-3) make list of stubs
  	std::unordered_map<int, std::vector<int>> stubs;
	for(int v=0; v<genG.N; ++v){
		int d = assigned_degree[v];
		if(stubs.find(d) == stubs.end()){
			stubs[d] = std::vector<int>();
		}
		int remain_stubs = assigned_degree[v] - int(genG.nlist[v].size());
		for(int i=0; i<remain_stubs; ++i){
			stubs[d].push_back(v);
		}
	}

	for(auto itr=stubs.begin(); itr!=stubs.end(); ++itr){
		int d = itr->first;
		std::shuffle(stubs[d].begin(), stubs[d].end(), engine);
	}

	//(4-4) connect each free stub of nodes with degree k and degree k' uniformly at random
	int u, v;
  	std::vector<int> ks;
  	target_num_jnt_deg.get_keys(ks);
  	std::shuffle(ks.begin(), ks.end(), engine);
  	for(int k:ks){
  		std::vector<int> ls;
  		target_num_jnt_deg.get_keys(ls);
  		std::shuffle(ls.begin(), ls.end(), engine);
  		for(int l:ls){
  			while(target_num_jnt_deg.get(k, l) != current_num_jnt_deg.get(k, l)){
				u = stubs[k].back();
		  		stubs[k].pop_back();
		  		v = stubs[l].back();
		  		stubs[l].pop_back();
		  		genG.add_edge(u,v);
		  		current_num_jnt_deg.add(k, l, 1);
		  		current_num_jnt_deg.add(l, k, 1);
			}
  		}
  	}

	genG.M = 0;
  	genG.maxd = 0;
  	for(i=0;i<genG.N;++i){
  		genG.M += genG.nlist[i].size();
  		if(int(genG.nlist[i].size()) > genG.maxd){
  			genG.maxd = genG.nlist[i].size();
  		}
  	}
  	genG.M = int(genG.M)/2;

	//(5) perform the targeting rewiring process so that a generated graph has approximately estimated degree-dependent clustering coefficient.
	//(5-1) estimate degree-dependent clustering coefficient
	Vector<double> est_ddcc;
	degree_dependent_clustering_coefficient_estimator(samplinglist, est_ddcc);

	std::vector<double> target_ddcc;
	int id, d;
	for(id=0; id<int(est_ddcc.entries.size()); ++id){
		d = est_ddcc.id_to_keys[id];
		if(d+1 > int(target_ddcc.size())){
			target_ddcc.resize(d+1, 0.0);
		}
		target_ddcc[d] = est_ddcc.get(d);
	}

	//(5-2) construct a set of rewirable edges
	std::vector<std::pair<int, int>> rewirable_edges;
	for(i=genG.n_q; i<int(genG.N); ++i){
		for(int j:genG.nlist[i]){
			if(j >= i && j >= genG.n_q){
				rewirable_edges.push_back({i, j});
			}
		}
	}

	//(5-3) targeting rewiring process
	Rewiring rew;
	rew.targeting_two_five_K(genG, target_ddcc, rewirable_edges);
	
  	return genG;
}