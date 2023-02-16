#include "cond_mutual_info.hpp"
#include <numeric>
#include <vector>

extern "C" {

info::dblvecvec *make_dblvecvec() { return new info::dblvecvec(); }

void attach_dblvec(
    info::dblvecvec *vec, const double vals[], const std::size_t length) {
	vec->emplace_back();
	auto &sub_vec = vec->back();
	sub_vec.reserve(length);

	for (std::size_t i = 0; i < length; ++i) {
		sub_vec.push_back(vals[i]);
	}
}

void delete_dblvecvec(info::cdblvecvec *vec) { delete vec; }

void cond_mutual_info(info::cdblvecvec *xs, info::cdblvecvec *ys,
    info::cdblvecvec *zs, const std::size_t p_samples, const double base,
    double *cmi_value, double *p_value) {
	auto cmi_estimator = info::cond_mutual_info(*xs, *ys, *zs);
	const auto result = cmi_estimator.calculate(p_samples, base);
	*cmi_value = result.first;
	*p_value = result.second;
}
}