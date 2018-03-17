/*************************************************************************
	> File Name: ARMAMath.h
	> Author: 
	> Mail: 
	> Created Time: 2018年03月15日 星期四 13时44分00秒
 ************************************************************************/

#ifndef _ARMAMATH_H
#define _ARMAMATH_H

#include <vector>

class ARMAMath{
public:
    double avgData(std::vector<double> dataArray);
    double sumData(std::vector<double> dataArray);
    double stderrData(std::vector<double> dataArray);
    double varerrData(std::vector<double> dataArray);
    std::vector<double>  autocorData(std::vector<double> dataArray,int order);
    std::vector<double>  autocovData(std::vector<double> dataArray,int order);

    double mutalCorr(std::vector<double> dataFir,std::vector<double> dataSec);
    double getModelAIC(std::vector<std::vector<double>> vec,std::vector<double> data,int type);

    std::vector<std::vector<double>> LevinsonSolve(std::vector<double> garma);
    std::vector<double> computeARCoe(std::vector<double> dataArray,int p);
    std::vector<double> computeMACoe(std::vector<double> dataArray,int q);
    std::vector<double> computeARMACoe(std::vector<double> dataArray,int p,int q);

};
#endif
