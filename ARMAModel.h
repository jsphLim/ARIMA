//
// Created by yue on 18-3-16.
//

#ifndef HUAWEI_ARMAMODEL_H
#define HUAWEI_ARMAMODEL_H


#include <vector>
#include "ARMAMath.cpp"

class ARMAModel{
private:
    std::vector<double> data;
    int p;
    int q;

public:
    ARMAModel(std::vector<double> data, int p,int q){
        this->data=data;
        this->p=p;
        this->q=q;
    }

    std::vector<std::vector<double>> solveCoeOfARMA(){
        std::vector<std::vector<double>> vec;
        ARMAMath ar_math;

        std::vector<double> armaCoe(ar_math.computeARMACoe(this->data,p,q));
        std::vector<double> arCoe(this->p+1);
        for(int i=0;i<arCoe.size();i++) arCoe[i]=armaCoe[i];

        std::vector<double>  maCoe(this->q+1);

        for(int i=0;i<maCoe.size();i++) {
            maCoe[i] = armaCoe[i+this->p+1];
        }
        vec.push_back(arCoe);
        vec.push_back(maCoe);
        //std::cout<<"fuck"<<std::endl;
       // std::cout<<vec.size()<<std::endl;
//        for(int i=0;i<vec.size();i++){
//            //for(int j=0;j<vec[i].size();j++)
//                 std::cout<<vec[i].size()<<std::endl;
//         }

        return vec;
    }
};
#endif //HUAWEI_ARMAMODEL_H
