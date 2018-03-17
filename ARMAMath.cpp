//
// Created by yue on 18-3-15.
//

#include <cstring>
#include <cmath>
#include <cstdlib>
#include "ARMAMath.h"

double ARMAMath::avgData(std::vector<double> dataArray) {
    return this->sumData(dataArray)/dataArray.size();
}


double ARMAMath::sumData(std::vector<double> dataArray) {
    double sumData=0;
    for(int i=0;i<dataArray.size();i++) sumData+=dataArray[i];

    return sumData;
}
double ARMAMath::stderrData(std::vector<double> dataArray) {
    return std::sqrt(this->varerrData(dataArray));
}


double ARMAMath::varerrData(std::vector<double> dataArray) {
    if(dataArray.size()<=1) return 0.0;
    double variance=0;
    double avgsumData=this->avgData(dataArray);

    for(int i=0;i<dataArray.size();i++){
        dataArray[i]-=avgsumData;
        variance+=dataArray[i]*dataArray[i];
    }
    return variance/(dataArray.size()-1);

}



std::vector<double>  ARMAMath::autocorData(std::vector<double> dataArray, int order) {
    std::vector<double>  autoCor;
    std::vector<double>  autoCov(this->autocovData(dataArray,order));
    double varData=this->varerrData(dataArray);
    if(varData!=0) {
        for (int i = 0; i < order ; i++) {
            autoCor[i]=autoCov[i]/varData;
        }
    }
    return autoCor;
}

std::vector<double>  ARMAMath::autocovData(std::vector<double> dataArray, int order) {
    std::vector<double>  autoCov(order+1);
    double mu = this->avgData(dataArray);
    for(int i=0;i<=order;i++){
        autoCov[i]=0.0;
        for(int j=0;j<dataArray.size()-i;j++){
            autoCov[i]+=(dataArray[j+i]-mu)*(dataArray[j]-mu);
        }
        autoCov[i]/=(dataArray.size()-i);
    }
    return autoCov;
}


double ARMAMath::mutalCorr(std::vector<double> dataFir, std::vector<double> dataSec) {
    double sumX=0.0;
    double sumY=0.0;
    double sumXY=0.0;
    double sumXSq=0.0;
    double sumYSq=0.0;
    int len=0;

    if(dataFir.size()!=dataSec.size()) len= (int) std::min(dataFir.size(), dataSec.size());
    else len= static_cast<int>(dataFir.size());

    for(int i=0;i<len;i++){
        sumX+=dataFir[i];
        sumY+=dataSec[i];
        sumXY+=dataFir[i]*dataSec[i];
        sumXSq+=dataFir[i]*dataFir[i];
        sumYSq+=dataSec[i]*dataSec[i];
    }

    double numerator = sumXY - sumX*sumY/len;
    double denominator = std::sqrt((sumXSq-sumX*sumX/len)*(sumYSq-sumY*sumY/len));

    if(denominator == 0) return 0.0;
    return numerator/denominator;

}

double gaussrand0()
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

double ARMAMath::getModelAIC(std::vector<std::vector<double>> vec, std::vector<double> data, int type) {
    int n = data.size();
    int p=0;
    int q=0;
    double tmpAR = 0.0,tmpMA=0.0;
    double sumErr=0.0;

    if(type==1){
        std::vector<double> maCoe = vec[0];
        q=(int)maCoe.size();
        std::vector<double> errData(q);

        for(int i=q-1;i<n;i++){
            tmpMA=0.0;
            for(int j=1;j<q;j++){
                tmpMA+=maCoe[j]*errData[j];
            }
            for(int j=q-1;j>0;j--){
                errData[j]=errData[j-1];
            }
            errData[0]=gaussrand0()*std::sqrt(maCoe[0]);
            sumErr+=(data[i]-tmpMA)*(data[i]-tmpMA);
        }
        return (n-(q-1))*std::log(sumErr/(n-(q-1)))+(q+1)*2;
    }
    else if(type==2){
        std::vector<double> arCoe = vec[0];
        p = (int)arCoe.size();

        for (int i = p - 1; i < n; ++i)
        {
            tmpAR = 0.0;
            for (int j = 0; j < p - 1; ++j)
            {
                tmpAR += arCoe[j] * data[i - j - 1];
            }
            sumErr += (data[i] - tmpAR) * (data[i] - tmpAR);
        }
//			return Math.log(sumErr) + (p + 1) * 2 / n;
        return (n - (p - 1)) * std::log(sumErr / (n - (p - 1))) + (p + 1) * 2;
    }
    else{
        std::vector<double> arCoe =vec[0];
        std::vector<double> maCoe=vec[1];
        p = (int)arCoe.size();
        q = (int)maCoe.size();
        std::vector<double> errData(q);

        for (int i = p - 1; i < n; ++i)
        {
            tmpAR = 0.0;
            for (int j = 0; j < p - 1; ++j)
            {
                tmpAR += arCoe[j] * data[i - j - 1];
            }
            tmpMA = 0.0;
            for (int j = 1; j < q; ++j)
            {
                tmpMA += maCoe[j] * errData[j];
            }

            for (int j = q - 1; j > 0; --j)
            {
                errData[j] = errData[j - 1];
            }
            errData[0] = gaussrand0() * std::sqrt(maCoe[0]);

            sumErr += (data[i] - tmpAR - tmpMA) * (data[i] - tmpAR - tmpMA);
        }
//			return Math.log(sumErr) + (q + p + 1) * 2 / n;
        return (n - (q + p - 1)) * std::log(sumErr / (n - (q + p - 1))) + (p + q) * 2;
    }
}

std::vector<std::vector<double>> ARMAMath::LevinsonSolve(std::vector<double> garma) {
    int order = garma.size()-1;
    std::vector<std::vector<double>> result;

    result.resize(order+1);
    for(int i=0;i<order+1;i++) result[i].resize(order+1);


    std::vector<double> sigmaSq(order+1);

    sigmaSq[0] = garma[0];
    result[1][1] = garma[1] / sigmaSq[0];
    sigmaSq[1] = sigmaSq[0] * (1.0 - result[1][1] * result[1][1]);
    for (int k = 1; k < order; ++k)
    {
        double sumTop = 0.0;
        double sumSub = 0.0;
        for (int j = 1; j <= k; ++j)
        {
            sumTop += garma[k + 1 - j] * result[k][j];
            sumSub += garma[j] * result[k][j];
        }
        result[k + 1][k + 1] = (garma[k + 1] - sumTop) / (garma[0] - sumSub);
        for (int j = 1; j <= k; ++j)
        {
            result[k + 1][j] = result[k][j] - result[k + 1][k + 1] * result[k][k + 1 - j];
        }
        sigmaSq[k + 1] = sigmaSq[k] * (1.0 - result[k + 1][k + 1] * result[k + 1][k + 1]);
    }
    result[0] = sigmaSq;

    return result;
}

std::vector<double> ARMAMath::computeARCoe(std::vector<double> dataArray, int p) {
    std::vector<double> garma = this->autocovData(dataArray,p);

    std::vector<std::vector<double>> result(this->LevinsonSolve(garma));

    std::vector<double> ARCoe(p+1);

    for(int i=0;i<p;i++){
        ARCoe[i] = result[p][i+1];

    }
    ARCoe[p] = result[0][p];
    return ARCoe;

}

std::vector<double> ARMAMath::computeMACoe(std::vector<double> dataArray, int q) {

    int p = (int)std::log(dataArray.size());

//		System.out.println("The best p is " + p);
    // 求取系数
    std::vector<double> bestGarma(this->autocovData(dataArray,p));
    std::vector<std::vector<double>> bestResult(this->LevinsonSolve(bestGarma));

    std::vector<double> alpha(p+1);
    alpha[0] = -1;
    for (int i = 1; i <= p; ++i)
    {
        alpha[i] = bestResult[p][i];
    }

    std::vector<double> paraGarma(q+1);
    for (int k = 0; k <= q; ++k)
    {
        double sum = 0.0;
        for (int j = 0; j <= p - k; ++j)
        {
            sum += alpha[j] * alpha[k + j];
        }
        paraGarma[k] = sum / bestResult[0][p];
    }

    std::vector<std::vector<double>> tmp (this->LevinsonSolve(paraGarma));
    std::vector<double> MACoe(q+1);
    for (int i = 1; i < MACoe.size(); ++i)
    {
        MACoe[i] = -tmp[q][i];
    }
    MACoe[0] = 1 / tmp[0][q];		//噪声参数

    return MACoe;
}

std::vector<double> ARMAMath::computeARMACoe(std::vector<double> dataArray, int p, int q) {
    std::vector<double> allGarma(this->autocovData(dataArray, p + q));
    std::vector<double> garma(p + 1);
    for (int i = 0; i < garma.size(); ++i)
    {
        garma[i] = allGarma[q + i];
    }
    std::vector<std::vector<double>> arResult(this->LevinsonSolve(garma));

    // AR
    std::vector<double> ARCoe(p+1);
    for (int i = 0; i < p; ++i)
    {
        ARCoe[i] = arResult[p][i + 1];
    }
    ARCoe[p] = arResult[0][p];
//		double [] ARCoe = this.YWSolve(garma);

    // MA
    std::vector<double> alpha(p+1);
    alpha[0] = -1;
    for (int i = 1; i <= p; ++i)
    {
        alpha[i] = ARCoe[i - 1];
    }

    std::vector<double> paraGarma(q+1);
    for (int k = 0; k <= q; ++k)
    {
        double sum = 0.0;
        for (int i = 0; i <= p; ++i)
        {
            for (int j = 0; j <= p; ++j)
            {
                sum += alpha[i] * alpha[j] * allGarma[std::abs(k + i - j)];
            }
        }
        paraGarma[k] = sum;
    }
    std::vector<std::vector<double>> maResult (this->LevinsonSolve(paraGarma));
    std::vector<double> MACoe(q+1);
    for (int i = 1; i <= q; ++i)
    {
        MACoe[i] = maResult[q][i];
    }
    MACoe[0] = maResult[0][q];

//		double [] tmp = this.YWSolve(paraGarma);
//		double [] MACoe = new double[q + 1];
//		System.arraycopy(tmp, 0, MACoe, 1, tmp.length - 1);
//		MACoe[0] = tmp[tmp.length - 1];

    std::vector<double> ARMACoe(p + q + 2);
    for (int i = 0; i < ARMACoe.size(); ++i)
    {
        if (i < ARCoe.size())
        {
            ARMACoe[i] = ARCoe[i];
        }
        else
        {
            ARMACoe[i] = MACoe[i - ARCoe.size()];
        }
    }
    return ARMACoe;
}
