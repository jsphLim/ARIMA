//
// Created by yue on 18-3-16.
//

#ifndef HUAWEI_ARIMAMODEL_H
#define HUAWEI_ARIMAMODEL_H

#include <vector>
#include "ARModel.h"
#include "MAModel.cpp"
#include "ARMAModel.h"

class ARIMAModel{
private:
    std::vector<double> dataArray;
    std::vector<double> dataFirDiff;

    std::vector<std::vector<double>> arima;
public:
    ARIMAModel(std::vector<double> dataArray) {this->dataArray.assign(dataArray.begin(),dataArray.end());}

    std::vector<double> preFirDiff(std::vector<double> preData){
        std::vector<double> res;
        for(int i=0;i<preData.size()-1;i++){
            double tmpData =preData[i+1]-preData[i];
            res.push_back(tmpData);
        }
        return res;
    }

    std::vector<double> preSeasonDiff(std::vector<double> preData){
        std::vector<double> res;

        for(int i=0;i<preData.size()-7;i++){

            double tmpData=preData[i+7]-preData[i];
            res.push_back(tmpData);
        }
        return res;
    }
    std::vector<double> preDealDiff(int period){
        if(period>=dataArray.size()-1){
            period=0;
        }

        switch (period){
            case 0: {
                return this->dataArray;
            }
            break;
            case 1: {
                std::vector<double> tmp(this->preFirDiff(this->dataFirDiff));
                this->dataFirDiff.assign(tmp.begin(), tmp.end());
                return this->dataFirDiff;
            }
            break;
            default: {
                return preSeasonDiff(dataArray);
            }
            break;
        }
    }

    std::vector<int> getARIMAModel(int period,std::vector<std::vector<int>> notModel,bool needNot){

        std::vector<double> data = this->preDealDiff(period);
        //for(int i=0;i<data.size();i++) std::cout<<data[i]<<std::endl;

        double minAIC = 1.7976931348623157E308D;
        std::vector<int> bestModel(3);
        int type = 0;
        std::vector<std::vector<double>> coe;


        // model产生, 即产生相应的p, q参数
        int len = data.size();

        if (len > 5)
        {
            len = 5;
        }




        int size = ((len + 2) * (len + 1)) / 2 - 1;
        std::vector<std::vector<int>> model;
        model.resize(size);
        for(int i=0;i<size;i++) model[i].resize(size);

        int cnt = 0;
        for (int i = 0; i <= len; ++i)
        {
            for (int j = 0; j <= len - i; ++j)
            {
                if (i == 0 && j == 0)
                    continue;
                model[cnt][0] = i;
                model[cnt++][1] = j;
            }
        }
        //std::cout<<size<<std::endl;
        for (int i = 0; i < cnt; ++i)
        {
            // 控制选择的参数

            bool token = false;
            if (needNot)
            {
                for (int k = 0; k < notModel.size(); ++k)
                {
                    if (model[i][0] == notModel[k][0] && model[i][1] == notModel[k][1])
                    {
                        token = true;
                        break;
                    }
                }
            }
            if (token)
            {
                continue;
            }

            if (model[i][0] == 0)
            {
                MAMoel* ma = new MAMoel(data,model[i][1]);

                //std::vector<std::vector<double>>
                coe=ma->solveCoeOfMA();
               // std::cout<<i<<coe.size()<<std::endl;
                //for(int ks=0;ks<ma->solveCoeOfMA().size();ks++) tmp.push_back(ma->solveCoeOfMA()[ks]);
                //coe.assign(tmp.begin(),tmp.end());
                type = 1;
            }
            else if (model[i][1] == 0)
            {
                ARModel* ar = new ARModel(data, model[i][0]);
                //std::vector<std::vector<double>> tmp(
                coe=ar->solveCoeOfAR();
             //   std::cout<<i<<coe.size()<<std::endl;
                //for(int ks=0;ks<ar->solveCoeOfAR().size();ks++) tmp.push_back(ar->solveCoeOfAR()[ks]);
                //coe.assign(tmp.begin(),tmp.end());
                type = 2;
            }
            else
            {
                //std::cout<<i<<model[i][0]<<" "<<model[i][1]<<std::endl;
                ARMAModel* arma = new ARMAModel(data, model[i][0], model[i][1]);;

                //std::vector<std::vector<double>> tmp(
                coe=arma->solveCoeOfARMA();
              //  std::cout<<i<<coe.size()<<std::endl;
                //for(int ks=0;ks<arma->solveCoeOfARMA().size();ks++) tmp.push_back(arma->solveCoeOfARMA()[ks]);
                //coe.assign(tmp.begin(),tmp.end());
                type = 3;
            }
            ARMAMath ar_math;
            double aic = ar_math.getModelAIC(coe, data, type);
            //std::cout<<aic<<std::endl;
            // 在求解过程中如果阶数选取过长，可能会出现NAN或者无穷大的情况

            if (aic<=1.7976931348623157E308D && !std::isnan(aic) && aic < minAIC)
            {
                minAIC = aic;
               // std::cout<<aic<<std::endl;
                bestModel[0] = model[i][0];
                bestModel[1] = model[i][1];
                bestModel[2] = (int)std::round(minAIC);
                this->arima = coe;
            }
        }
        return bestModel;
    }

    int aftDeal(int predictValue,int period){
        if (period >= dataArray.size())
        {
            period = 0;
        }

        switch (period)
        {
            case 0:
                return (int)predictValue;
            case 1:
                return (int)(predictValue + dataArray[dataArray.size() - 1]);
            case 2:
            default:
                return (int)(predictValue + dataArray[dataArray.size() - 7]);
        }
    }

    double gaussrand()
    {
        static double V1, V2, S;
        static int phase = 0;
        double X;

        if ( phase == 0 ) {
            do {
                double U1 = (double)rand() / RAND_MAX;
                double U2 = (double)rand() / RAND_MAX;

                V1 = 2 * U1 - 1;
                V2 = 2 * U2 - 1;
                S = V1 * V1 + V2 * V2;
            } while(S >= 1 || S == 0);

            X = V1 * sqrt(-2 * log(S) / S);
        } else
            X = V2 * sqrt(-2 * log(S) / S);

        phase = 1 - phase;

        return X;
    }


    int predictValue(int p,int q,int period){
        std::vector<double> data(this->preDealDiff(period));
        int n = data.size();
        int predict = 0;
        double tmpAR = 0.0, tmpMA = 0.0;
        std::vector<double> errData(q + 1);

        if (p == 0)
        {
            std::vector<double> maCoe(this->arima[0]);
            for(int k = q; k < n; ++k)
            {
                tmpMA = 0;
                for(int i = 1; i <= q; ++i)
                {
                    tmpMA += maCoe[i] * errData[i];
                }
                //产生各个时刻的噪声
                for(int j = q; j > 0; --j)
                {
                    errData[j] = errData[j - 1];
                }
                errData[0] = gaussrand()*std::sqrt(maCoe[0]);
            }

            predict = (int)(tmpMA); //产生预测
        }
        else if (q == 0)
        {
            std::vector<double> arCoe(this->arima[0]);

            for(int k = p; k < n; ++k)
            {
                tmpAR = 0;
                for(int i = 0; i < p; ++i)
                {
                    tmpAR += arCoe[i] * data[k - i - 1];
                }
            }
            predict = (int)(tmpAR);
        }
        else
        {
            std::vector<double> arCoe(this->arima[0]);
            std::vector<double> maCoe(this->arima[1]);

            for(int k = p; k < n; ++k)
            {
                tmpAR = 0;
                tmpMA = 0;
                for(int i = 0; i < p; ++i)
                {
                    tmpAR += arCoe[i] * data[k- i - 1];
                }
                for(int i = 1; i <= q; ++i)
                {
                    tmpMA += maCoe[i] * errData[i];
                }

                //产生各个时刻的噪声
                for(int j = q; j > 0; --j)
                {
                    errData[j] = errData[j-1];
                }

                errData[0] = gaussrand() * std::sqrt(maCoe[0]);
            }

            predict = (int)(tmpAR + tmpMA);
        }

        return predict;
    }

};
#endif //HUAWEI_ARIMAMODEL_H
