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
#include <iostream>
#include <iterator>

#include "double_utils.h"
#include "piecewise_function.h"
//#include <Rcpp.h>


/*
 * Printing coefficients
 */
std::ostream& operator<<(std::ostream &os, Coefficients const &c){
    return os << "(a=" << c.quadratic << ", b=" << c.linear << ", c=" << c.constant << ")";
}


/*
 * Function polling
 */
inline double gradient_at(Coefficients F, double x){
    if(greater(F.quadratic, 0)){
        return 2 * F.quadratic * x + F.linear;
    }
    else{
        return F.linear;
    }
}

inline bool is_increasing_at(Coefficients F, double x){
    return greater(gradient_at(F, x), 0);
}

inline bool is_decreasing_at(Coefficients F, double x){
    return less(gradient_at(F, x), 0);
}

inline double get_min(Coefficients F){
    if(!equal(F.quadratic, 0)){
        return -F.linear / (2 * F.quadratic);
    }
    else if(!equal(F.linear, 0)){
        return F.linear * -INFINITY;
    }
    else{
        // Flat function: minimum is everywhere
        return 1;
    }
}

inline bool min_in_interval(Coefficients F, double x1, double x2){
    // Note: the interval is ]x1, x2]
    if(equal(F.quadratic, 0) && equal(F.linear, 0)){
        return true;
    }
    else{
        double min = get_min(F);
        return less(x1, min) && (less(min, x2) || equal(x2, min));
    }
}


/*
 * Breakpoint tools
 */
inline double PiecewiseFunction::get_breakpoint_position(breakpoint_list_t::iterator b_ptr) {
    return (is_end(b_ptr)) ? INFINITY : b_ptr->first;
}

inline bool PiecewiseFunction::is_end(breakpoint_list_t::iterator b_ptr) {
    return b_ptr == this->breakpoint_coefficients.end();
}

inline bool PiecewiseFunction::is_begin(breakpoint_list_t::iterator b_ptr) {
    return b_ptr == this->breakpoint_coefficients.begin();
}


/*
 * The global minimum
 */
double PiecewiseFunction::get_minimum_position() {
    // Find the position of the minimum segment's minimum
    if(this->breakpoint_coefficients.size() == 0){
        return -INFINITY;
    }
    else if(equal(this->min_coefficients.quadratic, 0) && equal(this->min_coefficients.linear, 0)){
        if(is_end(this->min_ptr)){
            // Case: \___x lower bounds only
            return get_breakpoint_position(std::prev(this->min_ptr));
        }
        else if(is_begin(this->min_ptr)){
            // Case: ___x/ upper bounds only
            return get_breakpoint_position(this->min_ptr);
        }
        else{
            // Case: \__x__/
            return (get_breakpoint_position(std::prev(this->min_ptr)) + get_breakpoint_position(this->min_ptr)) / 2;
        }
    }
    else if(is_decreasing_at(this->min_coefficients, get_breakpoint_position(this->min_ptr))){
        // \x/ case
        return get_breakpoint_position(this->min_ptr);
    }
    else{
        // Case: there exists a single minimum and it is in the minimum segment
        return get_min(this->min_coefficients);
    }
}

double PiecewiseFunction::get_minimum_value() {
    double min_pos = this->get_minimum_position();
    if(std::isinf(min_pos)) {  // Only happens when there are no breakpoints in the function
        return 0;
    }
    return this->min_coefficients.quadratic * min_pos * min_pos + this->min_coefficients.linear * min_pos + this->min_coefficients.constant;
}


/*
 * Solver dynamic programming updates
 */
bool has_slack_at_position(double func_breakpoint_pos, double func_sign, Coefficients F, double pos){
    if(func_sign == -1 && less(pos, func_breakpoint_pos)){
        return true;
    }
    else if(func_sign == 1 && greater(pos, func_breakpoint_pos)){
        return true;
    }
    else{
        return false;
    }
}

int PiecewiseFunction::insert_point(double b, Coefficients F, bool is_upper_bound) {
    int n_pointer_moves = 0;

    double s;
    if(is_upper_bound){
        s = 1.;
    }
    else{
        s = -1.;
    }
    Coefficients Fs = F * s;

    if(this->breakpoint_coefficients.empty()){
        // Initialization
        breakpoint_list_t::iterator insert = this->breakpoint_coefficients.insert(breakpoint_t(b, Fs)).first;
        if(equal(s, 1)){
            this->min_ptr = insert;
            n_pointer_moves++;
        }
    }
    else{
        // General case

        // Insert the breakpoint
        std::pair<breakpoint_list_t::iterator, bool> insert = this->breakpoint_coefficients.insert(breakpoint_t(b, Fs));
        auto breakpoint_ptr = insert.first;
        auto breakpoint_was_added = insert.second;

        // If the breakpoint already exists, increase all its coefficients
        if (!breakpoint_was_added){
            breakpoint_ptr->second += Fs;
        }

        // Update the minimum function
        double min_pos = get_breakpoint_position(this->min_ptr);
        Coefficients G = this->min_coefficients;
        breakpoint_list_t::iterator g = this->min_ptr;
        if(has_slack_at_position(b, s, F, min_pos) || (equal(b, min_pos) && equal(s, -1))){
            G += F;
        }

        // Move the minimum pointer
        if(is_increasing_at(G, min_pos)){
            // Move left
            while(!is_begin(g) &&
                  !is_decreasing_at(G, min_pos) &&
                  !min_in_interval(G, get_breakpoint_position(std::prev(g)), min_pos)){
                g = std::prev(g);
                G -= g->second;
                min_pos = get_breakpoint_position(g);
                n_pointer_moves++;
            }
        }
        else if(is_decreasing_at(G, min_pos)){
            // Move right
            while(!is_end(g) &&
                  (min_in_interval(G + g->second, get_breakpoint_position(g), get_breakpoint_position(std::next(g)))  // If the next segment contains the minimum then move.
                   || !is_increasing_at(G + g->second, get_breakpoint_position(std::next(g))))){  // Otherwise if the next segment is not increasing, we'll eventually reach the minimum.
                G += g->second;
                g = std::next(g);
                min_pos = get_breakpoint_position(g);
                n_pointer_moves++;
            }
        }
        this->min_coefficients = G;
        this->min_ptr = g;
   }

    // Log progress
    if(this->verbose){
        /*std::cout << "\n\nINSERTION COMPLETED\n----------------------------" << std::endl;
      std::cout << "Breakpoint position: " << b << std::endl;
      std::cout << "Bound type: " << (is_upper_bound ? "upper" : "lower") << std::endl;
      std::cout << "N pointer moves: " << n_pointer_moves << std::endl;
      std::cout << "Minimum value: " << this->get_minimum_value() << std::endl;
      std::cout << "Minimum position: " << this->get_minimum_position() << std::endl;
      std::cout << "Minimum coefficients: " << this->min_coefficients << std::endl;
      std::cout << "The minimum pointer points to the breakpoint at " << get_breakpoint_position(this->min_ptr) << std::endl;
      std::cout << "Current breakpoints: {";
        for(auto bkpt : this->breakpoint_coefficients){
          std::cout << "(" << bkpt.first << ", " << bkpt.second << "), ";
        }
        std::cout << "}" << std::endl;
        std::cout << "\n\n" << std::endl;
        std::cout << "\n\n" << std::endl;*/
    }

    return n_pointer_moves;
}
