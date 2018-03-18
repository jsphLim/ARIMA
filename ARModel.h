//
// Created by yue on 18-3-16.
//

#ifndef ARMODEL_H
#define ARMODEL_H

#include <vector>
#include "ARMAMath.h"

class ARModel{
private:
    std::vector<double> data;
    int p;

public:
    ARModel(std::vector<double> data,int p){
        this->data.assign(data.begin(),data.end());
        this->p=p;
    }

    std::vector<std::vector<double> > solveCoeOfAR(){
        std::vector<std::vector<double> > vec;
        ARMAMath ar_math;
        std::vector<double>  arCoe(ar_math.computeARCoe(this->data,this->p));
        vec.push_back(arCoe);
        return vec;
    }
};


#endif //ARMODEL_H
