#ifndef CORE_DOUBLE_UTILS_H
#define CORE_DOUBLE_UTILS_H

#include <cmath>
#include <cstdlib>
#include <iostream>

#define TOL 1e-9

inline bool equal(double x, double y){
    if(std::isinf(x) || std::isinf(y)){
        return x == y;
    }
    else{
        return std::abs(x - y) <= TOL;
    }
}

inline bool not_equal(double x, double y){
    return !equal(x, y);
}

inline bool greater(double x, double y){
    if(std::isinf(x) || std::isinf(y)){
        return x > y;
    }
    else{
        return !equal(x, y) && x > y;
    }
}

inline bool less(double x, double y){
    if(std::isinf(x) || std::isinf(y)){
        return x < y;
    }
    else{
        return !equal(x, y) && x < y;
    }
}

class DoubleComparatorLess : public std::binary_function<double,double,bool>
{
public:
    bool operator()( const double &left, const double &right  ) const
    {
        return less(left, right);
    }
};

#endif //CORE_DOUBLE_UTILS_H
