//
// Created by yue on 18-3-16.
//

#ifndef HUAWEI_MAMODEL_H
#define HUAWEI_MAMODEL_H

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
#endif //HUAWEI_MAMODEL_H
