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
#ifndef MMIT_COEFFICIENTS_H
#define MMIT_COEFFICIENTS_H

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


#endif //MMIT_COEFFICIENTS_H
