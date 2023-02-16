#ifndef COND_MUTUAL_INFO_HPP
#define COND_MUTUAL_INFO_HPP
#include "common.hpp"
#include "discrete_dist.hpp"
#include "estimator.hpp"

namespace info {

/**
 * @brief A discrete conditional mutual information estimator: I(X;Y|Z).
 *
 */
class cond_mutual_info : public info::estimator {
	private:
	cdblvecvec &m_xs;
	cdblvecvec &m_ys;
	cdblvecvec &m_zs;
	info::discrete_dist m_xyz_dist;
	info::discrete_dist m_z_dist;
	info::discrete_dist m_xz_dist;
	info::discrete_dist m_yz_dist;

	/**
	 * @brief Reset the distributions with the data provided. When
	 * calculating p-values, shuffled versions of m_xs or m_ys may be passed
	 * into this function.
	 *
	 * @param xs
	 * @param ys
	 * @param zs
	 */
	inline void reset_dist(cdblvecvec &xs, cdblvecvec &ys, cdblvecvec &zs);
	inline void verify_data_integrity() const;

	/**
	 * @brief Calculate the conditional mutual information based on the
	 * current distributions
	 *
	 * @param base the log-base of the mutual information
	 * @return double
	 */
	double calculate_cmi(const double base) const;

	public:
	/**
	 * @brief Construct a new conditional discrete mutual information
	 * estimator: I(X;Y|Z).
	 *
	 * @param xs a vector of double vectors
	 * @param ys a vector of double vectors
	 * @param zs a vector of double vectors which the mutual information is
	 * conditional upon
	 */
	inline cond_mutual_info(cdblvecvec &xs, cdblvecvec &ys, cdblvecvec &zs);

	/**
	 * @brief Calculate the CMI and p-value
	 *
	 * @param p_samples the number of times the ys are shuffled to determine
	 * the p-value empirically
	 * @param base the base of the CMI value. The default is Euler's
	 * constant
	 * @return std::pair<double, double>
	 */
	std::pair<double, double> calculate(
	    const std::size_t p_samples, const double base) override;
};
} // namespace info

#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>

void info::cond_mutual_info::reset_dist(
    info::cdblvecvec &xs, info::cdblvecvec &ys, info::cdblvecvec &zs) {
	// assume data integrity was verified
	const std::size_t size = xs[0].size();
	const double sample_probability = (double)1.0 / size;
	const std::size_t var_count = xs.size() + ys.size() + zs.size();

	m_xyz_dist.m_probabilities.clear();
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

		m_xyz_dist.m_probabilities[event] += sample_probability;
	}

	std::vector<std::size_t> x_vars;
	std::vector<std::size_t> y_vars;
	std::vector<std::size_t> z_vars;

	// preallocate memory
	m_xyz_dist.m_variables.clear();
	m_xyz_dist.m_variables.reserve(var_count);
	x_vars.reserve(xs.size());
	y_vars.reserve(ys.size());
	z_vars.reserve(zs.size());

	for (std::size_t var_id = 0; var_id < var_count; ++var_id) {
		m_xyz_dist.m_variables.push_back(var_id);

		if (var_id < xs.size()) {
			x_vars.push_back(var_id);
		} else if (var_id < xs.size() + ys.size()) {
			y_vars.push_back(var_id);
		} else {
			z_vars.push_back(var_id);
		}
	}

	m_xyz_dist.reset_var_idx();

	m_xz_dist = m_xyz_dist.marginal(info::concat(x_vars, z_vars));
	m_yz_dist = m_xyz_dist.marginal(info::concat(y_vars, z_vars));
	m_z_dist = m_yz_dist.marginal(z_vars);
}

void info::cond_mutual_info::verify_data_integrity() const {
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

double info::cond_mutual_info::calculate_cmi(const double base) const {
	double cmi = 0.0;
	const double base_quotient = std::log(base);

	for (const auto &[xyz_e, xyz_p] : m_xyz_dist.m_probabilities) {

		const auto x_yz_e = info::split(xyz_e, m_xs.size());
		const auto y_z_e = info::split(x_yz_e.second, m_ys.size());

		const auto z_e = y_z_e.second;
		const auto xz_e = info::concat(x_yz_e.first, y_z_e.second);
		const auto yz_e = x_yz_e.second;

		const double z_p = m_z_dist.probability(z_e);
		const double xz_p = m_xz_dist.probability(xz_e);
		const double yz_p = m_yz_dist.probability(yz_e);

		cmi +=
		    xyz_p * std::log(z_p * xyz_p / xz_p / yz_p) / base_quotient;
	}

	return cmi;
}

info::cond_mutual_info::cond_mutual_info(
    info::cdblvecvec &xs, info::cdblvecvec &ys, info::cdblvecvec &zs)
    : m_xs(xs), m_ys(ys), m_zs(zs) {}

std::pair<double, double> info::cond_mutual_info::calculate(
    const std::size_t p_samples = 100,
    const double base = info::euler_constant) {
	verify_data_integrity();
	reset_dist(m_xs, m_ys, m_zs);
	const double cmi = calculate_cmi(base);

	if (p_samples == 0) {
		return std::make_pair(cmi, -1.0);
	}

	auto random_engine = std::default_random_engine{};
	auto ys_shuffle = m_ys;
	double p_val = 0.0;

	for (std::size_t rep = 0; rep < p_samples; ++rep) {

		for (auto &y : ys_shuffle) {
			std::shuffle(y.begin(), y.end(), random_engine);
		}

		reset_dist(m_xs, ys_shuffle, m_zs);

		if (calculate_cmi(base) >= cmi) {
			++p_val;
		}
	}

	p_val /= p_samples;

	return std::make_pair(cmi, p_val);
}

#endif // #ifndef COND_MUTUAL_INFO_HPP
