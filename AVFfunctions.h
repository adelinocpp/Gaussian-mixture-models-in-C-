#ifndef AVFFUNCTIONS_H
#define AVFFUNCTIONS_H

#include <string>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <complex>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <algorithm>
#include "matrix.h"

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

#include "gnuplot_i.hpp"

using namespace std;
using namespace math;

//--------------------------------------------------------------------------
struct MFCCparams
{
    bool bLogFile;          // Criar arquivo de log

    bool bComputeMFCC;      // Calcular as caracteristicas (MFCC) dos locutores

    bool bMultiThreading;   // Calcular utilizando várias threads

    int iMFCCsize;          // Numero de Componentes Mel-Cepstrais
    bool bVADselection;      // Seleciona os segmentos em que foram detectadas atividades de voz
    bool bFirstDD;          // Calcular primeira derivada temporal das MFCCs
    bool bSecondDD;         // Calcular segunda derivada temporal das MFCCs

    bool bUBM;              // Calcular Modelo de Fundo Universal
    bool bNormalizeByUBM;   // Normalizar Locutores pelo Modelo de Fundo Universal

    bool bNormalizeDD;      // Normalizar as derivadas pelo desvio padrao
    bool bPCA;              // Reduzir o MFCC pelos seus componentes principais.
    int iPCAsize;           // Numero de Componentes principais escolhidas.

    bool bComputeGMM;
    bool bComputeGMMofPCA;  // Computar GMM's dos componentes princiapais.
    int iMixSize;           // Tamanho do modelo de GMM
    bool bGMMbyGA;          // Calcular GMMs com algoritmo genético
    int iPopulationSize;    // Tamanho da populacao para algoritmo genetico
    bool bGMMbyEM;          // Calcular GMMs com maximizacao da expectancia
    int iEMmaxIteration;    // Numero maximo de iteracoes da maximizacao da expectancia
    bool bFuzzyKmeans;      // Utilizar clusterizacao Fuzzy K Means
    int iFKMmaxIteration;   // Numero maximo de iteracoes da clusterizacao Fuzzy K Means

    bool bComputeAllIndex;  // Calcular para todos sub-indices 1-Fala; 2-Palavras; 3-Frases; 4-Texto;
    int iComputeIndex;      // Especificar indice a ser calculado

    char * cAudioPath;      // Diretorio onde encontram-se os arquivos para processar
    char * cResultPath;     // Diretorio onde encontram-se os arquivos com os resultados
    char * cUBMPath;        // Diretorio onde encontra-se os arquivos do UBM
};
//--------------------------------------------------------------------------
struct ThreadParam
{
    int cpu;
    int busy;
    int * control;
    pthread_mutex_t * pMutex;

    string FileName;
    //string AudioPath;
    //string ResultPath;

    vector<double> v1;
    vector<double> v2;

    MFCCparams Parans;
};
//--------------------------------------------------------------------------------------------------

double minimun(matrix<double> mMatrix);
double maximun(matrix<double> mMatrix, int ModValue);
double minimun(vector<double> mMatrix);
double maximun(vector<double> mMatrix, int ModValue);
vector<double> minimun(matrix<double> mMatrix, int nDimension);
vector<double> maximun(matrix<double> mMatrix, int nDimension, int ModValue);

double maximun(matrix<double> mMatrix, int nDimension, int RowCol, int &Idx);

vector<int> IndexOfMax(matrix<double> mMatrix, int nDimension);
int IndexOfMax(vector<double> mVetor);


vector<double> mean(matrix<double> mMatrix, int nDimension);
vector<double> stdev(matrix<double> mMatrix, int nDimension);

matrix<double> correlation(matrix<double> mMatrix, int nDimension);
matrix<double> covariance(matrix<double> mMatrix, int nDimension);

matrix<double> StDNormalize(matrix<double> mMatrix, int nDimension, int afterIdx);

bool qr(matrix<double> A, matrix<double> &Q, matrix<double> &R);
bool svd(matrix<double> A, matrix<double> &U, matrix<double> &S);
bool pca(matrix<double> A, matrix<double> &PCA, int nComp);

void swap(double &a, double &b);
void swap(int &a, int &b);
int partition(vector<double> &vec, vector<int> &idx, int left, int right);
void quickSort(vector<double> &vec, vector<int> &idx, int left, int right);
void SortVector(vector<double> &vec, vector<int> &idx, int dec);
void SortByVector(vector<string> &vec, vector<int> idx, int dec);

void swap(string &a, string &b);
int partition(vector<string> &vec, vector<int> &idx, int left, int right);
void quickSort(vector<string> &vec, vector<int> &idx, int left, int right);
void SortVector(vector<string> &vec, vector<int> &idx, int dec);

vector<double> AcumFrequency(vector<double> mVector, int nBand);

vector<double> OverSample(vector<double> mVector, int nSamples, int nSideLobs);
vector<double> OverSample(vector<double> mVector, vector<double> &tVector, int nSamples, int nSideLobs);

vector<double> FindPeaks(vector<double> vec, vector<int> &idx, int dec);
vector<double> spline(vector<double> X, vector<double> Y, vector<double> Xx);

vector<double> UnwrapFunc(vector<double> vX, bool bRemTrend);

matrix<double> GetRandonSample(matrix<double> mX, int iSig, double fError);

matrix<double> csvread(string mFileName);
void csvwrite(matrix<double> mX, string mFileName);

double InformacaoMutua(matrix<double> mX);

//int NumSubElementos(matrix<double> mX, int iSec, matrix<double> &mY);
matrix<double> NumSubElementos(matrix<double> mX, int iSec, int &NumElementos);

// if dReadFileSize = true sort files by size; else sort by name
vector<string> read_directory(const string path, string sExt, bool dReadFileSize);

string RemoveExtension(string sNameFileExt);

double rand_normal(double mean, double stddev);


#endif // FUNCTIONS_H
