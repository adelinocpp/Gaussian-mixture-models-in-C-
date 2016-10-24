#ifndef CGMM_H
#define CGMM_H

#include "AVFfunctions.h"
//#include "CData.h"
//--------------------------------------------------------------------------
//--- Definicoes para HMM's
typedef matrix<double> mtxFeatures;
//--------------------------------------------------------------------------
//--- Definicoes paras GMM's
typedef matrix<double> mtxSig;
typedef vector<double> vecMu;
typedef vector<double> GeneVector;
//--------------------------------------------------------------------------
class CGMM
{
private:
    vector<double> mRho;
    vector<vecMu> mMu;
    vector<mtxSig>  mSigma;
    int nGaussians;
    bool CovDiag;
    bool Computed;
    unsigned long long int DataTag;

    //--- Variaveis para o Fuzzy EM
    int nMaxIterEM;
    int nMaxFKMIterEM;
    double dFKMFuziest;
    double epsFKM;
    double epsEM;
public:
    CGMM();
    ~CGMM() {};

    bool IsComputed(){return Computed;};
    int GetMixtureSize(){return mRho.size();};
    int GetMixtureDimension();
    vector<double> GetWeights(){return mRho;};
    vector<vecMu> GetCenters(){return mMu;};
    vector<mtxSig>  GetVariances(){return mSigma;};
    double GetWeight(int idx);
    vecMu GetCenter(int idx);
    mtxSig GetVariance(int idx);
    bool GetCovDiag(){return CovDiag;}
    unsigned long long int GetDataTag() {return DataTag;};
    int GetDimensionality(){return mMu[0].size();}

    void SetWeight(vector<double> vRho);
    void SetCenter(vector<vecMu> vMu);
    void SetVariance(vector<mtxSig>  vSigma);
    void SetWeight(double vRho, int idx);
    void SetCenter(vecMu vMu, int idx);
    void SetVariance(mtxSig vSigma, int idx);
    void SetCovDiag(bool bCovDiag){CovDiag = bCovDiag;};
    void SetMixtureSize(int nMix){nGaussians = nMix;};
    void SetComputed(bool bComputed){Computed = bComputed;};
    void SetDataTag(unsigned long long int ulliDataTag){DataTag = ulliDataTag;};
    void SetMaxEMIterations(int maxEM){nMaxIterEM = maxEM;};
    void SetMaxFKMIterations(int FKMIter){nMaxFKMIterEM = FKMIter;};

    //double Evaluate(vector<double> vX, int nMix,vector<double> extRho; vector<vecMu> extMu; vector<mtxSig>  extSigma);
    double EvaluateSumSigma(vector<double> vX, int nMix, mtxSig  extSigma);

    double Evaluate(double vX);
    double Evaluate(double vX, int i);
    double Evaluate(vector<double> vX);

    double Evaluate(vector<double> vX, int nMix);
    matrix<double> EvaluatePosteriori(matrix<double> vX);
    bool ComputeExpectationMaximization(matrix<double> vX, int nMixture, int nIter= 100, bool Fuzzy=true, int MaxFKIter=300);
    bool ComputeFuzzyKMeans(matrix<double> vX, int nMixture, double mFuzzy = 1.2, int nIter = 100, double eps = 1e-3);
    bool ImproveExpectationMaximization(matrix<double> vX, int nIter = 100, double eps = 1e-3);
    int GetMaxEMIterations(){return nMaxIterEM;};

    //--- Rotinas de leitura e escrita em arquivo
    vector<double> Encode();
    void Decode(vector<double> vGene, int nMixtures, int nDimensions);
    bool CSVRead(string mFileName);
    bool CSVWrite(string mFileName);

    double timeGMM;
};

//--------------------------------------------------------------------------
struct GMMparams
{
    bool bLogFile;

    int nGaussians;

    bool ComputeForEM;

    int nMaxIterEM;
    int nMaxFKMIterEM;

    int nIndividuals;
    int nGenerations;
    double pMutation;
};


#endif // CGMM_H
