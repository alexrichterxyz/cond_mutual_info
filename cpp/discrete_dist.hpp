#ifndef DISCRETE_DIST_HPP
#define DISCRETE_DIST_HPP
#include "common.hpp"
#include <unordered_map>
#include <vector>

namespace info {

class estimator;
class cond_mutual_info;

template <class T> struct vec_hash {
	std::hash<T> m_hash;
	inline std::size_t operator()(const std::vector<T> &vec) const noexcept;
};

using event_hash = vec_hash<double>;

class discrete_dist {
	private:
	// vector of variable ids (e.g. X_4, X_5 -> [4, 5])
	std::vector<std::size_t> m_variables;

	// mapping from variable id to index
	// the index determines the position in the event vector
	std::unordered_map<std::size_t, std::size_t> m_variables_idx;

	// mapping from event to probability (e.g. [X4=1.62, X_5=3.14] -> 0.5)
	std::unordered_map<dblvec, double, event_hash> m_probabilities;

	// create mapping from variable id to its index (i.e. reinitialize
	// m_variables_idx)
	inline void reset_var_idx();

	public:
	inline discrete_dist();

	/**
	 * @brief Create a distribution from double vectors. The index of the
	 * each double vector becomes the id of the corresponding random
	 * variable
	 *
	 * @param data, a vector of double vectors
	 */
	inline discrete_dist(cdblvecvec &data);

	/**
	 * @brief Get a reference to the variable ids (e.g. X_4, X_5 -> [4,
	 * 5])
	 *
	 * @return const std::vector<std::size_t>& reference to variable ids
	 */
	inline const std::vector<std::size_t> &variables() const;

	/**
	 * @brief Get a reference to the mapping from events to probabilites
	 * (e.g. [X4=1.62, X_5=3.14] -> 0.5)
	 *
	 * @return const std::unordered_map<dblvec, double,
	 * event_hash>& mapping from events to probabilities
	 */
	inline const std::unordered_map<dblvec, double, event_hash> &
	probabilities() const;

	/**
	 * @brief Get the probability of the specified event. Events have the
	 * form [(variable_id, value), ...]
	 *
	 * @param event, vector of pairs of the form [(variable_id, value), ...]
	 * @return double
	 */
	inline double probability(
	    const std::vector<std::pair<std::size_t, double>> &event);

	/**
	 * @brief Get the probability of an event
	 *
	 * @param event, a vector of values ordered by variable id. If the
	 * variable ids are [4, 5] and the event vector is [X_4=1.62, X_5=3.14],
	 * this corresponds to the event [X4=1.62, X_5=3.14]
	 * @return * double
	 */
	inline double probability(cdblvec &event) const;

	/**
	 * @brief Create a marginal distribution keeping only the variables with
	 * the ids specified
	 *
	 * @param keep, the vector of variable ids to keep
	 * @return discrete_dist
	 */
	inline discrete_dist marginal(const std::vector<std::size_t> &keep);

	/**
	 * @brief Create a conditional distribution. The condition has the form
	 * [(variable_id, value), ...]
	 *
	 * @param condition vector of pairs of the form [(variable_id, value),
	 * ...]
	 * @return discrete_dist
	 */
	inline discrete_dist conditional(
	    const std::vector<std::pair<std::size_t, double>> &condition);

	friend estimator;
	friend cond_mutual_info;
};

} // namespace info

#include <cassert>
#include <unordered_set>

template <class T>
std::size_t info::vec_hash<T>::operator()(
    const std::vector<T> &vec) const noexcept {
	std::size_t hash = 0x9e3779b9;

	for (const auto &v : vec) {
		// from boost
		hash ^= m_hash(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
	}

	return hash;
}

info::discrete_dist::discrete_dist() {}

info::discrete_dist::discrete_dist(info::cdblvecvec &data) {

	if (data.empty()) {
		return;
	}

	const std::size_t size = data[0].size();

	// check if all double vectors have the same size
	for (const auto &vec : data) {
		assert(vec.size() == size);
	}

	// assign distribution's variable ids
	m_variables.reserve(data.size());
	for (std::size_t var_id = 0; var_id < data.size(); ++var_id) {
		m_variables.push_back(var_id);
	}

	reset_var_idx();

	const double sample_probability = 1.0 / size;
	info::dblvec event(data.size(), 0.0); // reuse in each iteration
	for (std::size_t sample_idx = 0; sample_idx < size; ++sample_idx) {
		for (std::size_t var_idx = 0; var_idx < data.size();
		     ++var_idx) {
			event[var_idx] = data[var_idx][sample_idx];
		}

		m_probabilities[event] += sample_probability;
	}
}

void info::discrete_dist::reset_var_idx() {
	m_variables_idx.clear();
	m_variables.reserve(m_variables.size());

	for (std::size_t var_idx = 0; var_idx < m_variables.size(); ++var_idx) {
		m_variables_idx[m_variables[var_idx]] = var_idx;
	}
}

double info::discrete_dist::probability(
    const std::vector<std::pair<std::size_t, double>> &event) {
	double probability = 0.0;

	for (const auto &[obj_event, prob] : m_probabilities) {
		for (const auto &[var_id, val] : event) {
			if (obj_event[m_variables_idx[var_id]] != val) {
				goto event_dissatisfied;
			}
		}
		// event matched
		probability += prob;
	event_dissatisfied:;
	}

	return probability;
}

double info::discrete_dist::probability(info::cdblvec &event) const {
	const auto it = m_probabilities.find(event);

	if(it == m_probabilities.end()) {
		return 0.0;
	}

	return it->second;
}

const std::vector<std::size_t> &info::discrete_dist::variables() const {
	return m_variables;
}

const std::unordered_map<info::dblvec, double, info::event_hash> &
info::discrete_dist::probabilities() const {
	return m_probabilities;
}

info::discrete_dist info::discrete_dist::marginal(
    const std::vector<std::size_t> &keep) {
	discrete_dist new_dist;
	new_dist.m_variables = keep;
	new_dist.reset_var_idx();

	// determine marginal probabilities
	info::dblvec event(keep.size(), 0.0);
	for (const auto &[obj_event, probability] : m_probabilities) {

		std::size_t i = 0;
		for (const auto &var_id : keep) {
			event[i++] = obj_event[m_variables_idx[var_id]];
		}

		new_dist.m_probabilities[event] += probability;
	}

	return new_dist;
}

info::discrete_dist info::discrete_dist::conditional(
    const std::vector<std::pair<std::size_t, double>> &condition) {
	discrete_dist new_dist;

	// set variables; keep only those that are not in event
	new_dist.m_variables.reserve(m_variables.size() - condition.size());
	for (const auto &obj_var_id : m_variables) {
		for (const auto &[cond_var_id, _] : condition) {
			if (obj_var_id == cond_var_id) {
				goto do_not_keep;
			}
		}

		new_dist.m_variables.push_back(obj_var_id);
	do_not_keep:;
	}

	new_dist.reset_var_idx();

	// determine conditional probabilities
	double probability_sum = 0.0;
	dblvec event(new_dist.m_variables.size(), 0.0);
	for (const auto &[obj_event, probability] : m_probabilities) {
		// e = [2, 4, 1]
		std::size_t i = 0;

		// check if condition is satisfied
		for (const auto &[cond_var_id, cond_var_val] : condition) {
			if (obj_event[m_variables_idx[cond_var_id]] !=
			    cond_var_val) {
				goto condition_not_satisfied;
			}
		}

		// condition satisfied, create event
		for (const auto &var_id : new_dist.m_variables) {
			event[i++] = obj_event[m_variables_idx[var_id]];
		}

		probability_sum += probability;
		new_dist.m_probabilities[event] += probability;

	condition_not_satisfied:;
	}

	// normalize probabilities
	if (probability_sum > 0.0) {
		for (auto &[_, probability] : new_dist.m_probabilities) {
			probability /= probability_sum;
		}
	}

	return new_dist;
}

#endif // #ifndef DISCRETE_DIST_HPP
