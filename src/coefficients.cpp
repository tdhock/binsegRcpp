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
#include "coefficients.h"


Coefficients Coefficients::operator+(Coefficients &other) {
    return Coefficients(this->quadratic + other.quadratic, this->linear + other.linear, this->constant + other.constant);
}

void Coefficients::operator+=(Coefficients &other) {
    this->quadratic += other.quadratic;
    this->linear += other.linear;
    this->constant += other.constant;
}

Coefficients Coefficients::operator-(Coefficients &other) {
    return Coefficients(this->quadratic - other.quadratic, this->linear - other.linear, this->constant - other.constant);
}

void Coefficients::operator-=(Coefficients &other) {
    this->quadratic -= other.quadratic;
    this->linear -= other.linear;
    this->constant -= other.constant;
}

Coefficients Coefficients::operator*(double scalar) {
    return Coefficients(this->quadratic * scalar, this->linear * scalar, this->constant * scalar);
}

void Coefficients::operator*=(double scalar) {
    this->quadratic *= scalar;
    this->linear *= scalar;
    this->constant *= scalar;
}

Coefficients Coefficients::operator/(double scalar) {
    return Coefficients(this->quadratic / scalar, this->linear / scalar, this->constant / scalar);
}

void Coefficients::operator/=(double scalar) {
    this->quadratic /= scalar;
    this->linear /= scalar;
    this->constant /= scalar;
}

bool Coefficients::operator==(Coefficients &other){
    return this->quadratic == other.quadratic && this->linear == other.linear && this->constant == other.constant;
}