#include "conditional_mi.hpp"
#include <numeric>
#include <vector>

extern "C" {

std::vector<info::dblvec> *make_dblvec_vec() {
	return new std::vector<info::dblvec>();
}

void attach_dblvec(std::vector<info::dblvec> *vec, const double vals[],
    const std::size_t length) {
    vec->emplace_back();
    auto &sub_vec = vec->back();
    sub_vec.reserve(length);

    for (std::size_t i = 0; i < length; ++i) {
		sub_vec.push_back(vals[i]);
	}
}

void delete_dblvec_vec(const std::vector<info::dblvec> *vec) {
	delete vec;
}

void conditional_mi(const std::vector<info::dblvec> *xs, const std::vector<info::dblvec> *ys, const std::vector<info::dblvec> *zs, const std::size_t p_samples, double *cmi_value, double *p_value) {
    auto cmi_estimator = info::conditional_mi(*xs, *ys, *zs);
    const auto result = cmi_estimator.calculate(p_samples);
    *cmi_value = result.first;
    *p_value = result.second;
}

}