/*
MMIT: Max Margin Interval Trees
Copyright (C) 2017 Toby Dylan Hocking, Alexandre Drouin
This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#include <cmath>
#include <cstdlib>
#include <iostream>
#define TOL 1e-9
#include <map>
#include <set>

class Coefficients{
public:
    double quadratic;
    double linear;
    double constant;

    Coefficients(){this->quadratic = 0; this->linear = 0; this->constant = 0;}
    Coefficients(double a, double b, double c){this->quadratic = a; this->linear = b; this->constant = c;}
    Coefficients operator+(Coefficients &other);
    void operator+=(Coefficients &other);
    Coefficients operator-(Coefficients &other);
    void operator-=(Coefficients &other);
    Coefficients operator*(double scalar);
    void operator*=(double scalar);
    Coefficients operator/(double scalar);
    void operator/=(double scalar);
    bool operator==(Coefficients &other);
};

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

enum FunctionType{
    linear_hinge,
    squared_hinge
};

typedef std::map<double, Coefficients, DoubleComparatorLess> breakpoint_list_t;
typedef std::pair<double, Coefficients> breakpoint_t;

inline Coefficients get_l1_coefs
(const double y,
 const bool is_upper_limit){
  double s = (is_upper_limit ? 1. : -1.);
  return Coefficients(0, s, - s * y);
}

class PiecewiseFunction {

private:
    // Breakpoint and their coefficients
    breakpoint_list_t breakpoint_coefficients;

    // Minimum solution
    Coefficients min_coefficients;
    breakpoint_list_t::iterator min_ptr;  // Always on the right of the minimum

    // Utility vars + functions
    void construct(bool verbose){this->verbose = verbose; this->min_ptr = breakpoint_coefficients.end();}
    inline double get_breakpoint_position(breakpoint_list_t::iterator b_ptr);
    inline bool is_end(breakpoint_list_t::iterator b_ptr);
    inline bool is_begin(breakpoint_list_t::iterator b_ptr);
    bool verbose;

public:
    PiecewiseFunction(bool verbose){
        construct(verbose);
    }

    PiecewiseFunction(){
        construct(false);
    }

    // Point insertion
    int insert_point(double b, Coefficients F, bool is_upper_bound);
  void insert_l1(double data_value, double weight_value){
    insert_point
      (data_value,
       get_l1_coefs(data_value, false) * weight_value,
       false);
    insert_point
      (data_value,
       get_l1_coefs(data_value, true) * weight_value,
       true);
  }    

    // Minimum pointer functions
    double get_minimum_position();
    double get_minimum_value();
};



