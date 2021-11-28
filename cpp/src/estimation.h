#ifndef ESTIMATION_H
#define ESTIMATION_H

int size_estimator(const std::vector<SampledData> samplinglist, double &est_n);

int average_degree_estimator(const std::vector<SampledData> samplinglist, double &est_aved);

int degree_distribution_estimator(const std::vector<SampledData> samplinglist, Vector<double> &est_dd);

int jdd_estimator_ie(const std::vector<SampledData> samplinglist, const double est_n, const double est_aved, Matrix<double> &jdd_ie);

int jdd_estimator_te(const std::vector<SampledData> samplinglist, Matrix<double> &jdd_te);

int jdd_estimator_hybrid(const std::vector<SampledData> samplinglist, const double est_n, const double est_aved, Matrix<double> &jdd);

int degree_dependent_clustering_coefficient_estimator(const std::vector<SampledData> samplinglist, Vector<double> &est_ddcc);

#endif