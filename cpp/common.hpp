#ifndef INFO_COMMON_HPP
#define INFO_COMMON_HPP
#include <ostream>
#include <utility>
#include <vector>

namespace info {
using dblvec = std::vector<double>;
using cdblvec = const dblvec;

template <class T> T concat(const T &a, const T &b) {
	T ab;
	ab.reserve(a.size() + b.size());
	ab.insert(ab.end(), a.begin(), a.end());
	ab.insert(ab.end(), b.begin(), b.end());
	return ab;
}

// at argument specifies the index of the first element in the second object
template <class T> std::pair<T, T> split(const T &vec, const std::size_t at) {
	T a;
	T b;
	a.reserve(at);
	b.reserve(vec.size() - at);

	auto cutoff_it = vec.begin() + at;
	a.insert(a.end(), vec.begin(), cutoff_it);
	b.insert(b.end(), cutoff_it, vec.end());

	return std::make_pair(a, b);
}

std::ostream &operator<<(std::ostream &out, const std::vector<double> &vec) {
	out << "[ ";

	for (const auto &v : vec) {
		out << v << " ";
	}

	out << ']';
	return out;
}
} // namespace info

#endif // #ifndef INFO_COMMON_HPP