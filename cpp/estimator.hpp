#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP
#include <utility>

namespace info {

class estimator {
    public:
    virtual std::pair<double, double> calculate(std::size_t p_samples) = 0;
};

}

#endif // #ifndef ESTIMATOR_HPP