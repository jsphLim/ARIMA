//
// Created by yue on 18-3-16.
//

#ifndef MAMODEL_H
#define MAMODEL_H

#include <vector>
#include "ARMAMath.h"

class MAMoel{
private:
    std::vector<double> data;
    int p;

public:
    MAMoel(std::vector<double> data,int p);

    std::vector<std::vector<double>> solveCoeOfMA();
};
#endif //MAMODEL_H
