//
// Created by yue on 18-3-16.
//

#include "MAModel.h"


MAMoel::MAMoel(std::vector<double> data, int p) {
    this->data=data;
    this->p=p;
}

std::vector<std::vector<double>> MAMoel::solveCoeOfMA() {
    std::vector<std::vector<double> > vec;
    ARMAMath ar_math;
    std::vector<double>  maCoe(ar_math.computeMACoe(this->data,this->p));
    vec.push_back(maCoe);
    return vec;
}