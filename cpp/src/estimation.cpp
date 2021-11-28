#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "basic_function.h"
#include "graph.h"
#include "sampling.h"
#include "estimation.h"

int size_estimator(const std::vector<SampledData> samplinglist, double &est_n){
	int r = samplinglist.size();
	int k,l,v_k,v_l;
	int m = my_round(r*0.025); //parameter
	double d_k,d_l;
	double phi = 0.0; double psi = 0.0;

	for(k=0;k<r-m;++k){
		v_k = samplinglist[k].v;
		d_k = double(samplinglist[k].nlist.size());
		for(l=k+m;l<r;++l){
			v_l = samplinglist[l].v;
			d_l = double(samplinglist[l].nlist.size());
			if(v_k == v_l){
				phi += double(2);
			}
			psi += (double(d_k)/d_l)+ (double(d_l)/d_k);
		}
	}

	
	est_n = double(psi)/phi;
	
	return 0;
}

int average_degree_estimator(const std::vector<SampledData> samplinglist, double &est_aved){
	est_aved = 0.0;

	for(SampledData data:samplinglist){
		est_aved += double(1)/data.nlist.size();
	}

	est_aved = double(samplinglist.size())/est_aved;
	
	return 0;
}

int degree_distribution_estimator(const std::vector<SampledData> samplinglist, Vector<double> &est_dd){
	int i, d;
	double x = 0.0;
	est_dd = Vector<double>();

	for(SampledData data:samplinglist){
		d = data.nlist.size();
		est_dd.add(d, double(1)/d);
		x += double(1)/d;
	}

	for(i=0; i<int(est_dd.entries.size()); ++i){
		d = est_dd.id_to_keys[i];
		est_dd.set(d, double(est_dd.entries[i])/x);
	}

	return 0;
}

int jdd_estimator_ie(const std::vector<SampledData> samplinglist, const double est_n, const double est_aved, Matrix<double> &jdd_ie){
	//Induced Edges

	int r = samplinglist.size();
	int m = my_round(r*0.025);
	int i, j, v, k, l, c;
	double value;
	Matrix<double> phi;

	for(i=0;i<r-m;++i){
		v = samplinglist[i].v;
		k = samplinglist[i].nlist.size();
		for(j=i+m;j<r;++j){
			l = samplinglist[j].nlist.size();
			c = std::count(samplinglist[j].nlist.begin(), samplinglist[j].nlist.end(), v);
			value = double(c)/(k*l);
			phi.add(k, l, value);
			phi.add(l, k, value);
		}
	}

	jdd_ie = Matrix<double>();
	int num_sample = int((r-m)*(r-m+1));
	double sum_d = est_n*est_aved;
	for(i=0; i<int(phi.entries.size()); ++i){
		k = phi.id_to_keys[i];
		for(j=0; j<int(phi.entries[i].size()); ++j){
			l = phi.id_to_keys[j];
			value = double(phi.get(k, l))/num_sample;
			value *= sum_d;
			jdd_ie.set(k, l, value);
		}
	}

	return 0;
}

int jdd_estimator_te(const std::vector<SampledData> samplinglist, Matrix<double> &jdd_te){
	//Traversed Edges

	jdd_te = Matrix<double>();
	int r = samplinglist.size();
	int i, j, k,l;

	for(i=0;i<r-1;++i){
		k = samplinglist[i].nlist.size();
		l = samplinglist[i+1].nlist.size();
		jdd_te.add(k, l, 1);
		jdd_te.add(l, k, 1);
	}

	double value;
	for(i=0; i<int(jdd_te.entries.size()); ++i){
		k = jdd_te.id_to_keys[i];
		for(j=0; j<int(jdd_te.entries[i].size()); ++j){
			l = jdd_te.id_to_keys[j];
			value = double(jdd_te.get(k, l))/(2*(r-1));
			jdd_te.set(k, l, value);
		}
	}

	return 0;
}

int jdd_estimator_hybrid(const std::vector<SampledData> samplinglist, const double est_n, const double est_aved, Matrix<double> &jdd){
	Matrix<double> jdd_ie;
	Matrix<double> jdd_te;

	jdd_estimator_ie(samplinglist, est_n, est_aved, jdd_ie);
	jdd_estimator_te(samplinglist, jdd_te);

	jdd = Matrix<double>();

	int k, l;
	double value;
	for(int i=0; i<int(jdd_ie.entries.size()); ++i){
		k = jdd_ie.id_to_keys[i];
		for(int j=0; j<int(jdd_ie.entries[i].size()); ++j){
			l = jdd_ie.id_to_keys[j];
			if((k+l) >= 2*est_aved){
				value = jdd_ie.get(k, l);
				jdd.set(k, l, value);
				jdd.set(l, k, value);
			}
		}
	}

	for(int i=0; i<int(jdd_te.entries.size()); ++i){
		k = jdd_te.id_to_keys[i];
		for(int j=0; j<int(jdd_te.entries[i].size()); ++j){
			l = jdd_te.id_to_keys[j];
			if((k+l) < 2*est_aved){
				value = jdd_te.get(k, l);
				jdd.set(k, l, value);
				jdd.set(l, k, value);
			}
		}
	}

	return 0;
}

int degree_dependent_clustering_coefficient_estimator(const std::vector<SampledData> samplinglist, Vector<double> &est_ddcc){
	int r, i, c, d, v, s, t;
	r = samplinglist.size();
	Vector<double> phi;
	Vector<double> psi;

	for(i=0; i<r; ++i){
		d = samplinglist[i].nlist.size();
  		psi.add(d, double(1)/d);
	}

	for(i=1; i<r-1; ++i){
		d = samplinglist[i].nlist.size();
		if(d == 0 || d == 1){continue;}
		s = samplinglist[i-1].v;
		v = samplinglist[i].v;
		t = samplinglist[i+1].v;
		if(s != v && v != t && t != s){
			c = std::count(samplinglist[i+1].nlist.begin(), samplinglist[i+1].nlist.end(), s);
    		phi.add(d, double(c)/(d-1));
		}
	}

	for(i=0; i<int(phi.entries.size()); ++i){
		phi.entries[i] = double(phi.entries[i])/(r-2);
	}

	for(i=0; i<int(psi.entries.size()); ++i){
		psi.entries[i] = double(psi.entries[i])/r;
	}

	est_ddcc = Vector<double>();
	for(i=0; i<int(psi.entries.size()); ++i){
		d = psi.id_to_keys[i];
		if(psi.get(d) > 0){
			est_ddcc.set(d, double(phi.get(d))/psi.get(d));
		}
	}
	
	return 0;
}