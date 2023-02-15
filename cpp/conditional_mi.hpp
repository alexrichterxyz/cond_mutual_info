#ifndef CONDITIONAL_MI_HPP
#define CONDITIONAL_MI_HPP
#include "common.hpp"
#include "discrete_dist.hpp"
#include "estimator.hpp"

namespace info {

class conditional_mi : public info::estimator {
	private:
	const std::vector<dblvec> &m_xs;
	const std::vector<dblvec> &m_ys;
	const std::vector<dblvec> &m_zs;
	info::discrete_dist m_p_xyz;
	info::discrete_dist m_p_z;
	info::discrete_dist m_p_xz;
	info::discrete_dist m_p_yz;

	inline void reset_dist(const std::vector<dblvec> &xs,
    const std::vector<dblvec> &ys, const std::vector<dblvec> &zs);
	inline void verify_data_integrity() const;
	double calculate_cmi();
	
	public:
	inline conditional_mi(const std::vector<dblvec> &xs,
    const std::vector<dblvec> &ys, const std::vector<dblvec> &zs);
	std::pair<double, double> calculate(
	    const std::size_t p_samples) override;
};
} // namespace info

#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include <iostream>

void info::conditional_mi::reset_dist(const std::vector<dblvec> &xs,
    const std::vector<dblvec> &ys, const std::vector<dblvec> &zs) {
	const std::size_t size = xs[0].size();
	const double sample_probability = (double)1.0 / size;
	const std::size_t var_count = xs.size() + ys.size() + zs.size();

	m_p_xyz.m_probabilities.clear();
	info::dblvec event(var_count, 0.0);
	for (std::size_t sample_idx = 0; sample_idx < size; ++sample_idx) {
		std::size_t event_idx = 0;

		for (const auto &v : xs) {
			event[event_idx++] = v[sample_idx];
		}

		for (const auto &v : ys) {
			event[event_idx++] = v[sample_idx];
		}

		for (const auto &v : zs) {
			event[event_idx++] = v[sample_idx];
		}

		m_p_xyz.m_probabilities[event] += sample_probability;
	}

	m_p_xyz.m_variables.clear();
	std::vector<std::size_t> x_vars;
	std::vector<std::size_t> y_vars;
	std::vector<std::size_t> z_vars;

	// preallocate memory
	m_p_xyz.m_variables.clear();
	m_p_xyz.m_variables.reserve(var_count);
	x_vars.reserve(xs.size());
	y_vars.reserve(ys.size());
	z_vars.reserve(zs.size());

	for (std::size_t var_id = 0; var_id < var_count; ++var_id) {
		m_p_xyz.m_variables.push_back(var_id);

		if (var_id < xs.size()) {
			x_vars.push_back(var_id);
		} else if (var_id < xs.size() + ys.size()) {
			y_vars.push_back(var_id);
		} else {
			z_vars.push_back(var_id);
		}
	}

	m_p_xyz.reset_var_idx();

	m_p_xz = m_p_xyz.marginal(info::concat(x_vars, z_vars));
	m_p_yz = m_p_xyz.marginal(info::concat(y_vars, z_vars));
	m_p_z = m_p_xz.marginal(z_vars);
}

void info::conditional_mi::verify_data_integrity() const {
	assert(!m_xs.empty() && !m_ys.empty() && !m_zs.empty());

	// ensure every x, y, and z has same size
	const std::size_t size = m_xs[0].size();

	assert(size > 0);

	for (const auto &v : m_xs) {
		assert(v.size() == size);
	}

	for (const auto &v : m_ys) {
		assert(v.size() == size);
	}

	for (const auto &v : m_zs) {
		assert(v.size() == size);
	}
}

double info::conditional_mi::calculate_cmi() {
	double cmi = 0.0;
	
	for (const auto &[xyz_e, xyz_p] : m_p_xyz.m_probabilities) {

		const auto x_yz_e = info::split(xyz_e, m_xs.size());
		const auto y_z_e = info::split(x_yz_e.second, m_ys.size());

		const auto z_e = y_z_e.second;
		const auto xz_e = info::concat(x_yz_e.first, y_z_e.second);
		const auto yz_e = x_yz_e.second;

		const double z_p = m_p_z.probability(z_e);
		const double xz_p = m_p_xz.probability(xz_e);
		const double yz_p = m_p_yz.probability(yz_e);

		cmi += xyz_p * std::log(z_p * xyz_p / xz_p / yz_p);
	}

	return cmi;
}

info::conditional_mi::conditional_mi(const std::vector<dblvec> &xs,
    const std::vector<dblvec> &ys, const std::vector<dblvec> &zs)
    : m_xs(xs), m_ys(ys), m_zs(zs) {}

std::pair<double, double> info::conditional_mi::calculate(
    const std::size_t p_samples = 100) {
	verify_data_integrity();
	reset_dist(m_xs, m_ys, m_zs);
	const double cmi = calculate_cmi();
	
	if(p_samples == 0) {
		return std::make_pair(cmi, 0.0);
	}

	auto random_engine = std::default_random_engine {};
	auto ys_shuffle = m_ys;
	double p_val = 0.0;

	for(std::size_t rep = 0; rep < p_samples; ++rep) {
		
		for(auto &y: ys_shuffle){
			std::shuffle(y.begin(), y.end(), random_engine);
		}

		reset_dist(m_xs, ys_shuffle, m_zs);

		if(calculate_cmi() >= cmi) {
			++p_val;
		}
	}

	p_val /= p_samples;

	return std::make_pair(cmi, p_val);
}

#endif // #ifndef CONDITIONAL_MI_HPP
