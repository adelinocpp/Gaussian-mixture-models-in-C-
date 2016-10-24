#include "CGMM.h"

//--------------------------------------------------------------------------
CGMM::CGMM()
{
    mRho.clear();
    mMu.clear();
    mSigma.clear();
    nGaussians = 0;
    CovDiag = false; // Full covariace matrix by default
    Computed = false;
    DataTag = 0;
    timeGMM = 0;

    nMaxIterEM = 15;
    nMaxFKMIterEM = 50;
    dFKMFuziest = 1.1;
    epsFKM = 1e-3;
    epsEM = 1e-4;
}
//--------------------------------------------------------------------------
int CGMM::GetMixtureDimension()
{
    if (nGaussians == 0)
        return 0;
    else
        return mSigma[0].ColNo();
};
//--------------------------------------------------------------------------
double CGMM::GetWeight(int idx)
{
    if (idx < nGaussians)
        return mRho[idx];
    else
        return mRho[0];
};
//--------------------------------------------------------------------------
vecMu CGMM::GetCenter(int idx)
{
    if (idx < nGaussians)
        return mMu[idx];
    else
        return mMu[0];
};
//--------------------------------------------------------------------------
mtxSig CGMM::GetVariance(int idx)
{
    if (idx < nGaussians)
        return mSigma[idx];
    else
        return mSigma[0];
};
//--------------------------------------------------------------------------
void CGMM::SetWeight(vector<double> vRho)
{
    int i;
    nGaussians = vRho.size();
    mRho.clear();
    for(i = 0; i < nGaussians; i++)
        mRho.push_back(vRho[i]);
};
//--------------------------------------------------------------------------
void CGMM::SetCenter(vector<vecMu> vMu)
{
    int i;
    nGaussians = vMu.size();
    mMu.clear();
    for(i = 0; i < nGaussians; i++)
        mMu.push_back(vMu[i]);
};
//--------------------------------------------------------------------------
void CGMM::SetVariance(vector<mtxSig>  vSigma)
{
    int i;
    nGaussians = vSigma.size();
    mSigma.clear();
    for(i = 0; i < nGaussians; i++)
        mSigma.push_back(vSigma[i]);
};
//--------------------------------------------------------------------------
void CGMM::SetWeight(double vRho, int idx)
{
    if (idx > nGaussians)
        return;
    else
        mRho[idx] = vRho;
};
//--------------------------------------------------------------------------
void CGMM::SetCenter(vecMu vMu, int idx)
{
    if (idx > nGaussians)
        return;
    else
        mMu[idx] = vMu;
};
//--------------------------------------------------------------------------
void CGMM::SetVariance(mtxSig vSigma, int idx)
{
    if (idx > nGaussians)
        return;
    else
        mSigma[idx] = vSigma;
};
//--------------------------------------------------------------------------
double CGMM::EvaluateSumSigma(vector<double> vX, int nMix, mtxSig  extSigma)
{
    //  Avalia uma amostra de vX na gaussiana nMix da GMM soamndo ao Sigma original extSigma
    //  amostra dimensão 1 x n. n: dimensões.
    //  GMM n x 1. n: dimensões; e uma gaussiana.
    //  REtorna o valor da avaliação.
    int j, nDimensions;
    double SqrtDetSig, sumXmMUdSIG2, sumGauss;
    vecMu tMu;
    mtxSig tSig, invSig, XmMu, XmMuT, tRes;

    nDimensions = vX.size();
    if (nDimensions != this->GetDimensionality())
        return 0;
    /*
    if (nDimensions != extSigma.ColNo())
        return 0;
    if (mSigma.size() != extSigma.size())
        return 0;
    */

    if (nMix < 0)
        return 0;

    if (nMix < this->GetMixtureSize())
    {
        tMu = mMu[nMix];
        tSig = mSigma[nMix] + extSigma;
        sumXmMUdSIG2 = 0;
        SqrtDetSig = sqrt(tSig.Det());

        if (CovDiag)
            for(j = 0; j < nDimensions; j++)
                sumXmMUdSIG2 = sumXmMUdSIG2 + pow( (vX[j] - tMu[j]),2 )/(tSig(j,j));
        else
        {
            XmMuT.Null(1,nDimensions);
            XmMu.Null(nDimensions,1);;
            for(j = 0; j < 2; j++)
            {
                XmMuT(0,j) = vX[j] - tMu[j];
                XmMu(j,0) = vX[j] - tMu[j];
            }
            invSig = tSig.Inv();
            tRes = (XmMuT*invSig)*XmMu;
            sumXmMUdSIG2 = tRes(0,0);
            //sumXmMUdSIG2 = (double)((XmMuT*invSig)*XmMu);
            //printf("nxn: %i,%i, 0,0: %+5.2e\n",tRes(0,0),tRes.RowNo(), tRes.ColNo());
        }
        sumGauss = (mRho[nMix]/(pow(2*M_PI,((float)nDimensions/2))*SqrtDetSig))*exp(-0.5*sumXmMUdSIG2);
        return sumGauss;
    }
    else
    {
        sumGauss = 0;
        for (nMix = 0; nMix < this->GetMixtureSize();nMix++)
        {
            tMu = mMu[nMix];
            tSig = mSigma[nMix] + extSigma;
            sumXmMUdSIG2 = 0;
            SqrtDetSig = sqrt(tSig.Det());

            if (CovDiag)
                for(j = 0; j < nDimensions; j++)
                    sumXmMUdSIG2 = sumXmMUdSIG2 + pow( (vX[j] - tMu[j]),2 )/(tSig(j,j));
            else
            {
                XmMuT.Null(1,nDimensions);
                XmMu.Null(nDimensions,1);;
                for(j = 0; j < nDimensions; j++)
                {
                    XmMuT(0,j) = vX[j] - tMu[j];
                    XmMu(j,0) = vX[j] - tMu[j];
                }
                invSig = tSig.Inv();
                tRes = (XmMuT*invSig)*XmMu;
                sumXmMUdSIG2 = tRes(0,0);
                //printf("0,0: %+5.2e,%+5.2e  0,1 %+5.2e,%+5.2e 1,0: %+5.2e,%+5.2e 1,1: %+5.2e,%+5.2e\n",tRes(0,0),tRes(0,1),tRes(1,0),tRes(1,1));
            }
            sumGauss = sumGauss + (mRho[nMix]/(pow(2*M_PI,nDimensions/2)*SqrtDetSig))*exp(-0.5*sumXmMUdSIG2);
        }
        return sumGauss;
    }

    return sumGauss;

};
//--------------------------------------------------------------------------
double CGMM::Evaluate(double vX)
{
    //  Avalia uma amostra de vX no GMM (Modelo de mistura de gaussianas)
    //  amostra dimensão 1 x n. n: dimensões.
    //  GMM n x g. n: dimensões; e g: numero de gaussianas
    //  REtorna o valor da avaliação.
    int i, nMixtures;
    double sumXmMUdSIG2, sumGauss;
    vecMu tMu;
    mtxSig tSig;

    if (1 != this->GetDimensionality())
        return 0;

    nMixtures = this->GetMixtureSize();
    sumGauss = 0;
    for(i = 0; i < nMixtures; i++)
    {
        tMu.clear();
        tSig.Null();
        tMu = mMu[i];
        tSig = mSigma[i];

        sumXmMUdSIG2 = pow( (vX - tMu[0]),2)/(tSig(0,0));

        sumGauss = sumGauss + (mRho[i]/(sqrt(2*M_PI*tSig(0,0))))*exp(-0.5*sumXmMUdSIG2);
    }
    return sumGauss;
};
//--------------------------------------------------------------------------
double CGMM::Evaluate(double vX, int i)
{
    //  Avalia uma amostra de vX no GMM (Modelo de mistura de gaussianas)
    //  amostra dimensão 1 x n. n: dimensões.
    //  GMM n x g. n: dimensões; e g: numero de gaussianas
    //  REtorna o valor da avaliação.

    double sumGauss, sumXmMUdSIG2;
    vecMu tMu;
    mtxSig tSig;

    if (1 != this->GetDimensionality())
        return 0;

    tMu.clear();
    tSig.Null();
    tMu = mMu[i];
    tSig = mSigma[i];

    sumXmMUdSIG2 = pow( (vX - tMu[0]),2)/(tSig(0,0));

    sumGauss =  (mRho[i]/(sqrt(2*M_PI*tSig(0,0))))*exp(-0.5*sumXmMUdSIG2);

    return sumGauss;
};
//--------------------------------------------------------------------------
double CGMM::Evaluate(vector<double> vX)
{
    //  Avalia uma amostra de vX no GMM (Modelo de mistura de gaussianas)
    //  amostra dimensão 1 x n. n: dimensões.
    //  GMM n x g. n: dimensões; e g: numero de gaussianas
    //  REtorna o valor da avaliação.
    int i, j, nDimensions, nMixtures;
    double detSigma, sumXmMUdSIG2, sumGauss;
    vecMu tMu;
    mtxSig tSig, invSig, XmMu, XmMuT, tRes, DeltaSig;


    nDimensions = vX.size();
    if (nDimensions != this->GetDimensionality())
        return 0;

    nMixtures = this->GetMixtureSize();
    sumGauss = 0;
    for(i = 0; i < nMixtures; i++)
    {
        tMu.clear();
        tSig.Null();
        tMu = mMu[i];
        tSig = mSigma[i];
        detSigma = 1;
        sumXmMUdSIG2 = 0;
        if (CovDiag)
            for(j = 0; j < nDimensions; j++)
            {
                //sumXmMUdSIG2 = sumXmMUdSIG2 + pow( (vX[j] - tMu [j]),2)/(tSig(j,j)+1E-12);
                sumXmMUdSIG2 = sumXmMUdSIG2 + pow( (vX[j] - tMu [j]),2)/(tSig(j,j));
                detSigma = detSigma*tSig(j,j);
            }
        else
        {
            XmMuT.Null(1,nDimensions);
            XmMu.Null(nDimensions,1);;
            for(j = 0; j < nDimensions; j++)
            {
                XmMuT(0,j) = vX[j] - tMu[j];
                XmMu(j,0) = vX[j] - tMu[j];
            }

            invSig = tSig.Inv();

            tRes = (XmMuT*invSig)*XmMu;
            sumXmMUdSIG2 = tRes(0,0);
            detSigma = tSig.Det();
        }
        //detSigma = sqrt(abs(detSigma + 1E-9));
        detSigma = sqrt(abs(detSigma));
        sumGauss = sumGauss + (mRho[i]/(detSigma*pow(2*M_PI,nDimensions/2)))*exp(-0.5*sumXmMUdSIG2);
        if (isnan(sumGauss))
        {
            printf("Pressione para continuar...");
            getchar();
        }
    }
    return sumGauss;
};
//--------------------------------------------------------------------------
double CGMM::Evaluate(vector<double> vX, int nMix)
{
    //  Avalia uma amostra de vX na gaussiana nMix da GMM
    //  amostra dimensão 1 x n. n: dimensões.
    //  GMM n x 1. n: dimensões; e uma gaussiana.
    //  Retorna o valor da avaliação.
    int j, k, nDimensions, ic,il;
    double SqrtDetSig, sumXmMUdSIG2, sumGauss, Valor;
    vecMu tMu;
    mtxSig tSig, invSig, XmMu, XmMuT, tRes, DeltaSig;

    nDimensions = vX.size();
    if (nDimensions != this->GetDimensionality())
        return 0;

    if (nMix < 0)
        return 0;

    DeltaSig.Null(nDimensions,nDimensions);
    for (k = 0; k < nDimensions; k++)
        DeltaSig(k,k) = 0.1;

    if (nMix < this->GetMixtureSize())
    {
        tMu = mMu[nMix];
        tSig = mSigma[nMix];
        sumXmMUdSIG2 = 0;
        Valor = tSig.Det();
        if (Valor <= 0)
        {
            tSig = tSig + (tSig.Norm()/(nDimensions*nDimensions)) *DeltaSig;
            for (il = 0; il < nDimensions; il++)
            {
                for(ic = 0; ic < nDimensions; ic++)
                    printf("%+5.3f ",tSig(il,ic));
                printf("\n");
            }
            printf("Det: %+5.3f\n.",tSig.Det());
        }
        SqrtDetSig = sqrt(tSig.Det());
        //SqrtDetSig = sqrt(abs(tSig.Det()));
        invSig.Null(nDimensions,nDimensions);
        //printf("nxn: %i,%i, 0,0: %+5.2e\n",tSig.RowNo(), tSig.ColNo(), tSig(0,0));

        if (CovDiag)
            for(j = 0; j < nDimensions; j++)
                sumXmMUdSIG2 = sumXmMUdSIG2 + pow( (vX[j] - tMu[j]),2 )/(tSig(j,j));
        else
        {
            //printf("46 - %i; %i\n",nDimensions,vX.size());
            XmMuT.Null(1,nDimensions);
            XmMu.Null(nDimensions,1);;
            for(j = 0; j < nDimensions; j++)
            {
                XmMuT(0,j) = vX[j] - tMu[j];
                XmMu(j,0) = vX[j] - tMu[j];
            }

            invSig = tSig.Inv();
            tRes = (XmMuT*invSig);
            tRes = tRes*XmMu;
            sumXmMUdSIG2 = tRes(0,0);

            //sumXmMUdSIG2 = (double)((XmMuT*invSig)*XmMu);
            //printf("nxn: %i,%i, 0,0: %+5.2e\n",tRes(0,0),tRes.RowNo(), tRes.ColNo());
        }
        sumGauss = (mRho[nMix]/(pow(2*M_PI,((float)nDimensions/2))*SqrtDetSig))*exp(-0.5*sumXmMUdSIG2);

        if (std::isnan(sumGauss))
                printf("NaN in evaluate nMix: %i\n",nMix);
        return sumGauss;
    }
    else
    {
        sumGauss = 0;

        for (nMix = 0; nMix < this->GetMixtureSize();nMix++)
        {
            tMu = mMu[nMix];
            tSig = mSigma[nMix];
            sumXmMUdSIG2 = 0;
            SqrtDetSig = sqrt(tSig.Det());

            if (CovDiag)
                for(j = 0; j < nDimensions; j++)
                    sumXmMUdSIG2 = sumXmMUdSIG2 + pow( (vX[j] - tMu[j]),2 )/(tSig(j,j));
            else
            {
                XmMuT.Null(1,nDimensions);
                XmMu.Null(nDimensions,1);;
                for(j = 0; j < nDimensions; j++)
                {
                    XmMuT(0,j) = vX[j] - tMu[j];
                    XmMu(j,0) = vX[j] - tMu[j];
                }
                invSig = tSig.Inv();
                tRes = (XmMuT*invSig)*XmMu;
                sumXmMUdSIG2 = tRes(0,0);
                //printf("0,0: %+5.2e,%+5.2e  0,1 %+5.2e,%+5.2e 1,0: %+5.2e,%+5.2e 1,1: %+5.2e,%+5.2e\n",tRes(0,0),tRes(0,1),tRes(1,0),tRes(1,1));
                //sumXmMUdSIG2 = (double)((XmMuT*invSig)*XmMu);
            }
            sumGauss = sumGauss + (mRho[nMix]/(pow(2*M_PI,nDimensions/2)*SqrtDetSig))*exp(-0.5*sumXmMUdSIG2);
        }
       // if (std::isnan(sumGauss))
       //         printf("NaN in evaluate nMix: %i\n",nMix);
        return sumGauss;
    }

    return sumGauss;
};
//--------------------------------------------------------------------------
matrix<double> CGMM::EvaluatePosteriori(matrix<double> vX)
{
    //  Avalia as amostras de vX na GMM
    //  Amostra dimensão d x n. d: Numero de amostras; e n: dimensões.
    //  GMM n x g. n: dimensões; e g: o numero de gaussiana.
    //  Retorna o matriz g x d. g: numero de gaussianas; e d: o numero de amostras.

    int i, j, k, nDimensions, nSamples;
    vector<double> sX;
    matrix<double> pX;
    double sumPriori, vPriori;

    nSamples = vX.ColNo();
    nDimensions = vX.RowNo();

    pX.Null(nGaussians,nSamples);

    for(i = 0; i < nSamples; i++)
    {
        sX.clear();
        for (j = 0; j < nDimensions; j++)
        {
            //valor = vX(j,i);
            sX.push_back(vX(j,i));

        }


        sumPriori = 0;
        for (k = 0; k < nGaussians; k++)
        {
            vPriori = this->Evaluate(sX,k);
            pX(k,i) = vPriori;
            sumPriori = sumPriori + vPriori;
        }
        for (k = 0; k < nGaussians; k++)
            if (!(pX(k,i) == 0))
                pX(k,i) = pX(k,i)/(sumPriori);
    }
    return pX;
};
//--------------------------------------------------------------------------
bool CGMM::ComputeFuzzyKMeans(matrix<double> vX, int nMixture, double mFuzzy, int nIter, double eps)
{
    int xSize, xDimensions;
    int i, j, k, d, idx, lastCol;
    double SumUkRow, argRand, sumDen, sumNum, auxSum, fMaxDiff;
    matrix<double> Uk, Ukm, deltaUk, TempSigma;
    mtxFeatures tempX;
    vector<double> auxMu, auxNum, meanMu;
    vector<int> idxMaxUk;
    vector<mtxFeatures> vecSamplesByCluster;

    nGaussians = nMixture;
    xSize = vX.ColNo();
    xDimensions = vX.RowNo();

    //--- Iniciate Uk -------
    Uk.Null(nMixture,xSize);
    for(i = 0; i < xSize; i++)
    {
        SumUkRow = 0;
        for(j = 0; j < nMixture; j++)
        {
            argRand = ((double)rand()/RAND_MAX);
            Uk(j,i) = argRand;
            //Uk(j,i) = 1;
            SumUkRow = SumUkRow + Uk(j,i);
        }
        for(j = 0; j < nMixture; j++)
            Uk(j,i) = Uk(j,i)/SumUkRow;
    };

    Ukm.Null(nMixture,xSize);
    deltaUk.Null(nMixture,xSize);

    /*
    //---- rotina verifica Uk
    string sLine;
    char lBuff[64];
    float sumMix;
    for(i = 0; i < floor(xSize/20); i++)
    {
        sLine.clear();
        sumMix = 0;
        for(j = 0; j < nMixture; j++)
        {
            sprintf(lBuff,"%5.3f, ",Uk(j,i));
            sLine.append(lBuff);
            sumMix += Uk(j,i);
        }
        printf("%s = %5.3f\n",sLine.c_str(),sumMix);
    };
    printf("Verifica fuzzy k means encerrada, pressione ENTER para continuar... \n");
    std::cin.get();
    //---- fimrotina verifica Uk
    */

    idx = 0;
    while (idx < nIter)
    {
        mMu.clear();
        for(j = 0; j < nMixture; j++)
        {
            auxMu.clear();
            auxNum.clear();
            for(d = 0; d < xDimensions; d++)
            {
                auxNum.push_back(0.0);
                auxMu.push_back(0.0);
            }
            sumDen = 0;
            for(i = 0; i < xSize; i++)
            {
                sumDen = sumDen + pow(Uk(j,i),mFuzzy);
                for(d = 0; d < xDimensions; d++)
                {
                    auxNum[d] = auxNum[d] + pow(Uk(j,i),mFuzzy)*vX(d,i);
                }

            }

            for(d = 0; d < xDimensions; d++)
                auxMu[d] = auxNum[d]/sumDen;

            mMu.push_back(auxMu);
        }

        idx++;
        for(j = 0; j < nMixture; j++)
        {

            for(i = 0; i < xSize; i++)
            {
                auxSum = 0;
                for(d = 0; d < xDimensions; d++)
                    auxSum = auxSum + pow((mMu[j][d] - vX(d,i)),2);

                sumNum = sqrt(auxSum);  // "sumNum" is the distance between outMu[j] - vX(:,i)

                sumDen = 0;
                for(k = 0; k < nMixture; k++)
                {
                    auxSum = 0;
                    for(d = 0; d < xDimensions; d++)
                        auxSum = auxSum + pow((mMu[k][d] - vX(d,i)),2);

                    sumDen = sumDen + pow(sumNum/sqrt(auxSum),(2/(mFuzzy - 1)));
                }
                Ukm(j,i) = 1/sumDen;
            }

        }
        // --- Convergence test

        for(i = 0; i < xSize; i++)
            for(j = 0; j < nMixture; j++)
                deltaUk(j,i) = abs(Uk(j,i) - Ukm(j,i));

        fMaxDiff = maximun(deltaUk,0);
        if (fMaxDiff < eps)
            break;

        for(i = 0; i < xSize; i++)
            for(j = 0; j < nMixture; j++)
                Uk(j,i) = Ukm(j,i);
    }

    // Gera Rho e Sigma
    mRho.clear();
    mSigma.clear();

    idxMaxUk = IndexOfMax(Uk, 1);

    //---- fim rotina verifica Uk

    for (i = 0; i < nMixture; i++)
    {
        tempX.Null();
        tempX.SetSize(xDimensions,1);
        vecSamplesByCluster.push_back(tempX);
    }


    for (i = 0; i < xSize; i++)
    {
        idx = idxMaxUk[i];
        tempX = vecSamplesByCluster[idx];


        lastCol = tempX.ColNo() - 1;
        for (j = 0; j < xDimensions; j++)
                tempX(j,lastCol) = vX(j,i);

        tempX.SetSize(tempX.RowNo(),tempX.ColNo() + 1);
        vecSamplesByCluster[idx] = tempX;
    }


    for(idx = 0; idx < nMixture; idx++)
    {
        tempX = vecSamplesByCluster[idx];
        tempX.SetSize(tempX.RowNo(),tempX.ColNo() - 1);
        vecSamplesByCluster[idx] = tempX;
        d = tempX.ColNo();

        mRho.push_back((float)d/xSize);
        TempSigma.Null(xDimensions,xDimensions);
        TempSigma.SetSize(xDimensions,xDimensions);

        if (d > 3)
            TempSigma = covariance(tempX, 1);
        else
        {
            for(i = 0; i < xDimensions; i++)
                TempSigma(i,i) = 1;
        }
        mSigma.push_back(TempSigma);
        if (TempSigma.Det() <= 0)
        {
            int il, ic;
            printf("Determinante: %+5.4f, x(%lu,%lu) \n.",TempSigma.Det(),tempX.RowNo(),tempX.ColNo());
            for (il = 0; il < xDimensions; il++)
            {
                for(ic = 0; ic < xDimensions; ic++)
                    printf("%+5.3f ",TempSigma(il,ic));
                printf("\n");
            }
            printf("fim.\n");
        }
    }

    Uk.~matrix();
    Ukm.~matrix();
    return true;
};
//--------------------------------------------------------------------------
bool CGMM::ComputeExpectationMaximization(matrix<double> vX, int nMixture, int nIter, bool Fuzzy, int MaxFKIter)
{
    int j, k;
    int xDimensions;
    double tRho, sumRho, tRand, dSigma;
    vector<double> tempMu;
    matrix<double> tempSigma;
    vector<double> vXmax, vXmin, vXline, vXfreq;
    matrix<double> posterioriX;
    vector<double> sumPosX, sumPosX2; //Posteriori prob sum's

    xDimensions = vX.RowNo();
    vXmax = maximun(vX,1,0);
    vXmin = minimun(vX,1);

    mRho.clear();
    mMu.clear();
    mSigma.clear();
    nGaussians = nMixture;

    srand ( time(NULL) );

     //--- Inicializa as misturas
    if(Fuzzy)
    {
        ComputeFuzzyKMeans(vX, nMixture,dFKMFuziest,nIter,MaxFKIter);
    }
    else
    {
        sumRho = 0;
        for(k = 0; k < nMixture; k++)
        {
            tRho = ((double)rand()/RAND_MAX);
            sumRho = sumRho + tRho;
            mRho.push_back(tRho);
            tempMu.clear();
            tempSigma.Null(xDimensions,xDimensions);
            for(j = 0; j < xDimensions; j++)
            {
                tRand = ((double)rand()/RAND_MAX);
                tempMu.push_back( (0.1 + 0.9*tRand)*(vXmax[j]-vXmin[j]) + vXmin[j]);
                tRand = ((double)rand()/RAND_MAX);
                dSigma =  pow(abs(vXmax[j]-vXmin[j])/(6*nGaussians),2);
                tempSigma(j,j) = dSigma;
            }
            mMu.push_back(tempMu);
            mSigma.push_back(tempSigma);
        };
        for(k = 0; k < nMixture; k++)
            mRho[k] = mRho[k]/sumRho;
    }

    this->Computed = this->ImproveExpectationMaximization(vX,nIter);

    return Computed;
};
//--------------------------------------------------------------------------
bool CGMM::ImproveExpectationMaximization(matrix<double> vX, int nIter, double eps)
{
    int i, j, k, idxIter;
    int xSize, xDimensions, nMixture;
    vector<double> tempMu;
    matrix<double> tempSigma;
    matrix<double> posterioriX;
    double sumPosP, MaxDiffMu, MaxDiffSig, TMaxDiffMu, TMaxDiffSig, rLog, rLogAnt, sumT;
    vector<double> sumPosMu, sumPosDSig, diffMu; //Posteriori prob sum's
    matrix<double> sumPosFSig, XmMu, XmMuT, diffSig;

    vector<vecMu> AntMu;
    vector<mtxSig>  AntSigma;

    xSize = vX.ColNo();
    xDimensions = vX.RowNo();

    nMixture = mRho.size();
    nGaussians = nMixture;

    MaxDiffMu = 0;
    MaxDiffSig = 0;
    AntMu = mMu;
    AntSigma = mSigma;

    //tNCol = vX.ColNo();
    //tNRow = vX.RowNo();

    posterioriX = EvaluatePosteriori(vX);
    rLogAnt = 0;
    for(i = 0; i < xSize; i++)
    {
        sumT = 0;
        for(k = 0; k < nMixture; k++)
            sumT += posterioriX(k,i);
        if (sumT > 0)
            rLogAnt += log10(sumT);
    }

    for(idxIter = 0; idxIter < nIter; idxIter++)
    {
        posterioriX = EvaluatePosteriori(vX);

        mRho.clear();
        mMu.clear();
        mSigma.clear();
        rLog = 0;
        // start of maximization
        for(k = 0; k < nMixture; k++)
        {
            sumPosP = 0;
            sumPosMu.clear();
            sumPosDSig.clear();

            //---- Inicializa as variaveis
            if (CovDiag)    // (CovDiag == true) => Covariace matrix is diagonal
                for(j = 0; j < xDimensions; j++)
                {
                    sumPosMu.push_back(0.0);
                    sumPosDSig.push_back(0.0);
                }
            else            // (CovDiag == false) => Covariace matrix is full
            {
                sumPosFSig.Null(xDimensions,xDimensions);
                for(j = 0; j < xDimensions; j++)
                    sumPosMu.push_back(0.0);
            }

            //--- Computa parametros do modelo

            tempMu.clear();
            tempSigma.Null(xDimensions,xDimensions);
            if (CovDiag)    // (CovDiag == true) => Covariace matrix is diagonal
            {
                for(i = 0; i < xSize; i++)
                {
                    sumPosP = sumPosP + posterioriX(k,i);
                    for(j = 0; j < xDimensions; j++)
                        sumPosMu[j] = sumPosMu[j] + (double)posterioriX(k,i)*vX(j,i);
                }
                if (sumPosP == 0)
                    sumPosP += 1e-5;

                for(j = 0; j < xDimensions; j++)
                    tempMu.push_back(sumPosMu[j]/sumPosP);

                for(i = 0; i < xSize; i++)
                {
                    for(j = 0; j < xDimensions; j++)
                        sumPosDSig[j] = sumPosDSig[j] + (double)posterioriX(k,i)*(vX(j,i) - tempMu[j])*(vX(j,i) - tempMu[j]);
                }
                for(j = 0; j < xDimensions; j++)
                   tempSigma(j,j) = sumPosDSig[j]/sumPosP;
            }
            else            // (CovDiag == false) => Covariace matrix is full
            {
                for(i = 0; i < xSize; i++)
                {
                    sumPosP = sumPosP + posterioriX(k,i);
                    for(j = 0; j < xDimensions; j++)
                        sumPosMu[j] = sumPosMu[j] + (double)posterioriX(k,i)*vX(j,i);
                }

                if (sumPosP == 0)
                    sumPosP += 1e-5;

                for(j = 0; j < xDimensions; j++)
                    tempMu.push_back(sumPosMu[j]/sumPosP);

                sumPosFSig.Null(xDimensions,xDimensions);
                for(i = 0; i < xSize; i++)
                {
                    XmMu.Null(xDimensions,1);
                    XmMuT.Null(1,xDimensions);
                    for(j = 0; j < xDimensions; j++)
                    {
                        XmMu(j,0) = vX(j,i) - tempMu[j];
                        XmMuT(0,j) = vX(j,i) - tempMu[j];
                    }
                    sumPosFSig = sumPosFSig + posterioriX(k,i)*(XmMu*XmMuT);
                }
                tempSigma = sumPosFSig/sumPosP;
            }

            mRho.push_back(sumPosP/xSize);
            mMu.push_back(tempMu);
            mSigma.push_back(tempSigma);
        }
        // end of maximization


        // Convergence criteria
        MaxDiffMu = 0;
        MaxDiffSig = 0;
        for(k = 0; k < nMixture; k++)
        {
            diffMu.clear();
            for(j = 0; j < xDimensions; j++)
                diffMu.push_back(abs(AntMu[k][j] - mMu[k][j]));
            diffSig = AntSigma[k] - mSigma[k];
            TMaxDiffMu = maximun(diffMu,1);
            TMaxDiffSig = maximun(diffSig,1);
            if (TMaxDiffMu > MaxDiffMu)
                MaxDiffMu = TMaxDiffMu;
            if (TMaxDiffSig > MaxDiffSig)
                MaxDiffSig = TMaxDiffSig;
        }

        if ((MaxDiffMu < eps) && (MaxDiffSig < eps))
            break;

        posterioriX = EvaluatePosteriori(vX);
        for(i = 0; i < xSize; i++)
        {
            sumT = 0;
            for(k = 0; k < nMixture; k++)
                sumT += posterioriX(k,i);
            if (sumT > 0)
                rLog += log10(sumT);
        }

        AntMu = mMu;
        AntSigma = mSigma;
        rLogAnt = rLog;
        if (rLog < rLogAnt)
            ComputeFuzzyKMeans(vX, nMixture,dFKMFuziest,nMaxFKMIterEM,epsFKM);
    }
    return true;
};
//--------------------------------------------------------------------------
vector<double> CGMM::Encode()
{
    int i, j, k, nMixtures, nDimensions;
    vector<double> TempGene;

    TempGene.clear();
    nMixtures = this->GetMixtureSize();
    nDimensions = this->GetMixtureDimension();
    for(i = 0; i < nMixtures; i++)
        TempGene.push_back(this->mRho[i]);
    for(i = 0; i < nMixtures; i++)
        for(j = 0; j < nDimensions; j++)
            TempGene.push_back(this->mMu[i][j]);

    for(i = 0; i < nMixtures; i++)
        if (CovDiag)    // (CovDiag == true) => Covariace matrix is diagonal
        {
            for(j = 0; j < nDimensions; j++)
                TempGene.push_back(this->mSigma[i](j,j));
        }
        else            // (CovDiag == false) => Covariace matrix is full
        {
            for(j = 0; j < nDimensions; j++)
                for(k = j; k < nDimensions; k++)
                TempGene.push_back(this->mSigma[i](j,k));
        }

    return TempGene;
}
//--------------------------------------------------------------------------
void CGMM::Decode(vector<double> vGene, int nMixtures, int nDimensions)
{
    int i, j, k, idx;
    mtxSig mVariance;
    vecMu vCenter;
    this->nGaussians = nMixtures;

    mRho.clear();
    mMu.clear();
    mSigma.clear();
    for(i = 0; i < nMixtures; i++)
    {
        idx = i;
        this->mRho.push_back(vGene[i]);
    }

    for(i = 0; i < nMixtures; i++)
    {
        vCenter.clear();
        for(j = 0; j < nDimensions; j++)
        {
            idx = nMixtures + i*nDimensions + j;
            vCenter.push_back(vGene[idx]);
        }
        this->mMu.push_back(vCenter);
    }

    for(i = 0; i < nMixtures; i++)
    {
        mVariance.Null(nDimensions,nDimensions);
        if (CovDiag)    // (CovDiag == true) => Covariace matrix is diagonal
        {
            for(j = 0; j < nDimensions; j++)
            {
                idx = nMixtures*(nDimensions + 1) + i*nDimensions + j;
                mVariance(j,j) = vGene[idx];
            }
        }
        else            // (CovDiag == false) => Covariace matrix is full
        {
            for(j = 0; j < nDimensions; j++)
                for(k = j; k < nDimensions; k++)
                {
                    idx = nMixtures*(nDimensions + 1) + i*(int)((float)nDimensions*((float)nDimensions+1)/2) + j*nDimensions + k - j*(j+1)/2;     // Valor do indice????????????????
                    if (j == k)
                        mVariance(j,k) = vGene[idx];
                    else
                    {
                        mVariance(j,k) = vGene[idx];
                        mVariance(k,j) = vGene[idx];
                    }
                }
        }
        this->mSigma.push_back(mVariance);
    }
};
//--------------------------------------------------------------------------
bool CGMM::CSVRead(string mFileName)
{
    FILE * pFile;
    int i, nSamples, nMixtures, nDimensions, nCovDiag;
    double fTemp;
    vector<double> TempGene;

    pFile = fopen (mFileName.c_str(),"r");

    fscanf(pFile, "%llu",&DataTag);
    fscanf(pFile, "%lf",&timeGMM);
    fscanf(pFile, "%i",&nMixtures);
    fscanf(pFile, "%i",&nDimensions);
    fscanf(pFile, "%i",&nCovDiag);

    nGaussians = nMixtures;
   // printf("%i\n",nGaussians);
   // std::cin.get();
    if (nCovDiag == 1)  // (CovDiag == true) => Covariace matrix is diagonal
    {
        CovDiag = true;
        nSamples = nMixtures*(1 + 2*nDimensions);
    }
    else                // (CovDiag == false) => Covariace matrix is full
    {
        CovDiag = false;
        nSamples = nMixtures*(1 + nDimensions + nDimensions*(nDimensions+1)/2);
    }

    TempGene.clear();
    for (i = 0; i < nSamples; i++)
    {
        fscanf(pFile,"%lf",&fTemp);
        TempGene.push_back(fTemp);
    }
    this->Decode(TempGene, nMixtures, nDimensions);

    fclose (pFile);
    return true;
};
//--------------------------------------------------------------------------
bool CGMM::CSVWrite(string mFileName)
{
    FILE * pFile;
    int i, nSamples, nMixtures, nDimensions;
    string strTemp;
    char strValor[16];
    vector<double> TempGene;

    pFile = fopen(mFileName.c_str(),"w");

    nMixtures = this->GetMixtureSize();
    nDimensions = this->GetMixtureDimension();

    fprintf (pFile, "%llu\n",DataTag);
    fprintf (pFile, "%8.4e\n",timeGMM);
    fprintf (pFile, "%i\n",nMixtures);
    fprintf (pFile, "%i\n",nDimensions);

    if (CovDiag)    // (CovDiag == true) => Covariace matrix is diagonal
        fprintf (pFile, "1\n");
    else            // (CovDiag == false) => Covariace matrix is full
        fprintf (pFile, "0\n");

    TempGene.clear();
    TempGene = this->Encode();
    nSamples = TempGene.size();

    for (i = 0; i < nSamples; i++)
    {
        strTemp.clear();
        sprintf(strValor,"%6.5e ",TempGene[i]);
        strTemp.append(strValor);
        fprintf(pFile,"%s\n",strTemp.c_str());
    }

    fclose (pFile);
    return true;
};
//--------------------------------------------------------------------------
//  Fim CGMM
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
