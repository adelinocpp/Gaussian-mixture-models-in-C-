#include "AVFfunctions.h"

using namespace std;
using namespace math;
//--------------------------------------------------------------------------
double minimun(matrix<double> mMatrix)
{
    unsigned int i, j;
    double tMin;
    tMin = mMatrix(0,0);
    for(i = 0; i < mMatrix.ColNo(); i++)
        for(j = 0; j < mMatrix.RowNo(); j++)
            if (mMatrix(j,i) < tMin)
                tMin = mMatrix(j,i);

    return tMin;
};
//--------------------------------------------------------------------------
double maximun(matrix<double> mMatrix, int ModValue)
{
    unsigned int i, j;
    double tMax;
    tMax = mMatrix(0,0);

    if (ModValue == 0) // Compute normal values
    {
    for(i = 0; i < mMatrix.ColNo(); i++)
        for(j = 0; j < mMatrix.RowNo(); j++)
            if (mMatrix(j,i) > tMax)
                tMax = mMatrix(j,i);
    }
    else // Compute in modular values
    {
        for(i = 0; i < mMatrix.ColNo(); i++)
            for(j = 0; j < mMatrix.RowNo(); j++)
                if (abs(mMatrix(j,i)) > tMax)
                    tMax = abs(mMatrix(j,i));
    }
    return tMax;
};
//--------------------------------------------------------------------------
double minimun(vector<double> mMatrix)
{
    unsigned int i;
    double tMin;
    tMin = mMatrix[0];
    for(i = 0; i < mMatrix.size(); i++)
        if (mMatrix[i] < tMin)
            tMin = mMatrix[i];

    return tMin;
};
//--------------------------------------------------------------------------
double maximun(vector<double> mMatrix, int ModValue)
{
    unsigned int i;
    double tMax;
    tMax = mMatrix[0];
    if (ModValue == 0) // Compute normal values
    {
        for(i = 0; i < mMatrix.size(); i++)
            if (mMatrix[i] > tMax)
                tMax = mMatrix[i];
    }
    else // Compute in modular values
    for(i = 0; i < mMatrix.size(); i++)
            if (abs(mMatrix[i]) > tMax)
                tMax = abs(mMatrix[i]);
    return tMax;
};
//--------------------------------------------------------------------------
vector<double> minimun(matrix<double> mMatrix, int nDimension)
{
    // If nDimension = 0 compute de minimun of each column
    // else compute the minimun of each line
    unsigned int i, j;
    double tMin;
    vector<double> vMin;

    if (nDimension == 0)
    {
        for(i = 0; i < mMatrix.ColNo(); i++)
        {
            tMin = DBL_MAX;
            for(j = 0; j < mMatrix.RowNo(); j++)
                if (mMatrix(j,i) < tMin)
                    tMin = mMatrix(j,i);
            vMin.push_back(tMin);
        }
    }
    else
    {
        for(j = 0; j < mMatrix.RowNo(); j++)
        {
            tMin = DBL_MAX;
            for(i = 0; i < mMatrix.ColNo(); i++)
                if (mMatrix(j,i) < tMin)
                    tMin = mMatrix(j,i);

            vMin.push_back(tMin);
        }
    }
    return vMin;
};
//--------------------------------------------------------------------------
vector<double> maximun(matrix<double> mMatrix, int nDimension, int ModValue)
{
    // If nDimension = 0 compute de maximun of each column
    // else compute the maximun of each line
    unsigned int i, j;
    double tMax;
    vector<double> vMax;

    if (ModValue)
    {
        if (nDimension == 0)
        {
            for(i = 0; i < mMatrix.ColNo(); i++)
            {
                tMax = -DBL_MAX;
                for(j = 0; j < mMatrix.RowNo(); j++)
                    if (mMatrix(j,i) > tMax)
                        tMax = mMatrix(j,i);
                vMax.push_back(tMax);
            }
        }
        else
        {
            for(j = 0; j < mMatrix.RowNo(); j++)
            {
                tMax = -DBL_MAX;
                for(i = 0; i < mMatrix.ColNo(); i++)
                    if (mMatrix(j,i) > tMax)
                        tMax = mMatrix(j,i);
                vMax.push_back(tMax);
            }
        }
    }
    else
    {
        if (nDimension == 0)
        {
            for(i = 0; i < mMatrix.ColNo(); i++)
            {
                tMax = -DBL_MAX;
                for(j = 0; j < mMatrix.RowNo(); j++)
                    if (abs(mMatrix(j,i)) > tMax)
                        tMax = abs(mMatrix(j,i));
                vMax.push_back(tMax);
            }
        }
        else
        {
            for(j = 0; j < mMatrix.RowNo(); j++)
            {
                tMax = -DBL_MAX;
                for(i = 0; i < mMatrix.ColNo(); i++)
                    if (abs(mMatrix(j,i)) > tMax)
                        tMax = abs(mMatrix(j,i));
                vMax.push_back(tMax);
            }
        }
    }
    return vMax;
};
//---------------------------------------------------------------------
double maximun(matrix<double> mMatrix, int nDimension, int RowCol, int &Idx)
{
    // If nDimension = 0 compute de maximun of a row indicate by RowCol
    // else compute de maximun of a col indicate by RowCol
    unsigned int i;
    double tMax;
    double vMax;

    if (nDimension == 0)
    {
        Idx = -1;
        tMax = -DBL_MAX;
        for(i = 0; i < mMatrix.ColNo(); i++)
        {
            if (mMatrix(RowCol,i) > tMax)
            {
                    tMax = mMatrix(RowCol,i);
                    Idx = i;
            }

        }
        vMax = tMax;
    }
    else
    {
        Idx = -1;
        tMax = -DBL_MAX;
        for(i = 0; i < mMatrix.RowNo(); i++)
        {
            if (mMatrix(i,RowCol) > tMax)
            {
                tMax = mMatrix(i,RowCol);
                Idx = i;
            }
        }
        vMax = tMax;
    }
    return vMax;

};
//---------------------------------------------------------------------
vector<int> IndexOfMax(matrix<double> mMatrix, int nDimension)
{
    // -- If (nDimension = 0) compute the maximum of each row
    // return a vector of size RowNo
    // -- else compute the maximum of each col
    // return a vector of size ColNo
    unsigned int i, j, idxMax;
    double tMax;
    vector<int> vMax;

    vMax.clear();
    if (nDimension == 0)
    {
        for(j = 0; j < mMatrix.RowNo(); j++)
        {
            idxMax = -1;
            tMax = -DBL_MAX;
            for(i = 0; i < mMatrix.ColNo(); i++)
                if (mMatrix(j,i) > tMax)
                {
                    tMax = mMatrix(j,i);
                    idxMax = i;
                }
            vMax.push_back(idxMax);
        }
    }
    else
    {
        for(i = 0; i < mMatrix.ColNo(); i++)
        {
            idxMax = -1;
            tMax = -DBL_MAX;
            for(j = 0; j < mMatrix.RowNo(); j++)
                if (mMatrix(j,i) > tMax)
                {
                    tMax = mMatrix(j,i);
                    idxMax = j;
                }
            vMax.push_back(idxMax);
        }
    }
    return vMax;
};
//---------------------------------------------------------------------
int IndexOfMax(vector<double> mVetor)
{
    unsigned int i, idx_max;
    double tMax;
    tMax = 0;

    for(i = 0; i < mVetor.size(); i++)
        if (mVetor[i] > tMax)
            idx_max = i;

    return idx_max;
};
//---------------------------------------------------------------------
vector<double> mean(matrix<double> mMatrix, int nDimension)
{
    // If nDimension = 0 compute de mean of each column
    // else compute the maximun of each line
    unsigned int i, j;
    double tMed;
    vector<double> vMed;

    vMed.clear();
    if (nDimension == 0)
    {
        for(i = 0; i < mMatrix.ColNo(); i++)
        {
            tMed = 0;
            for(j = 0; j < mMatrix.RowNo(); j++)
                tMed += mMatrix(j,i);
            vMed.push_back(tMed/mMatrix.RowNo());
        }
    }
    else
    {
        for(j = 0; j < mMatrix.RowNo(); j++)
        {
            tMed = 0;
            for(i = 0; i < mMatrix.ColNo(); i++)
                tMed += mMatrix(j,i);
            vMed.push_back(tMed/mMatrix.ColNo());
        }
    }
    return vMed;
};
//---------------------------------------------------------------------
vector<double> stdev(matrix<double> mMatrix, int nDimension)
{
    // If nDimension = 0 compute de standard deviation of each column
    // else compute the maximun of each line
    unsigned int i, j;
    double tMed, tStdev;
    vector<double> vMed, vStdev;

    vMed.clear();
    vStdev.clear();
    if (nDimension == 0)
    {
        for(i = 0; i < mMatrix.ColNo(); i++)
        {
            tMed = 0;
            for(j = 0; j < mMatrix.RowNo(); j++)
                tMed += mMatrix(j,i);
            vMed.push_back(tMed/mMatrix.RowNo());
        }
        for(i = 0; i < mMatrix.ColNo(); i++)
        {
            tStdev = 0;
            for(j = 0; j < mMatrix.RowNo(); j++)
                tStdev += pow(mMatrix(j,i) - vMed[i],2);
            vStdev.push_back(sqrt(tStdev/(mMatrix.RowNo() -1)));
        }
    }
    else
    {
        for(j = 0; j < mMatrix.RowNo(); j++)
        {
            tMed = 0;
            for(i = 0; i < mMatrix.ColNo(); i++)
                tMed += mMatrix(j,i);
            vMed.push_back(tMed/mMatrix.ColNo());
        }
        for(j = 0; j < mMatrix.RowNo(); j++)
        {
            tStdev = 0;
            for(i = 0; i < mMatrix.ColNo(); i++)
                tMed += pow(mMatrix(j,i) - vMed[j],2);;
            vStdev.push_back(sqrt(tStdev/(mMatrix.ColNo() -1)));
        }
    }
    return vStdev;
};
//---------------------------------------------------------------------
matrix<double> correlation(matrix<double> mMatrix, int nDimension)
{
    // nDimension == 0; correlation by lines
    // nDimension <> 0; correlation by cols
    int i, j, k, nLin, nCol;
    vector<double> vMean, vStd;
    double tMean, tStd, xyStd;
    nLin = mMatrix.RowNo();
    nCol = mMatrix.ColNo();
    matrix<double> rMatrix;

    if(nDimension == 0)
    {
        rMatrix.Null(nCol,nCol);

        vMean.clear();
        for (i = 0; i < nCol; i++)
        {
            tMean = 0;
            for (j = 0; j < nLin; j++)
                tMean += mMatrix(j,i);
            vMean.push_back(tMean/nLin);
        }
        vStd.clear();
        for (i = 0; i < nCol; i++)
        {
            tStd = 0;
            for (j = 0; j < nLin; j++)
                tStd += pow( (mMatrix(j,i) - vMean[j]),2);
            vStd.push_back(sqrt(tStd));
        }
        for (i = 0; i < nCol; i++)
        {
            for (j = 0; j < nCol; j++)
            {
                xyStd = 0;
                for(k = 0; k < nLin; k++)
                    xyStd += (mMatrix(k,i) - vMean[i])*(mMatrix(k,j) - vMean[j]);
                rMatrix(i,j) = xyStd/(vStd[i]*vStd[j]);
            }

        }

    }
    else
    {
        rMatrix.Null(nLin,nLin);

        vMean.clear();
        for (i = 0; i < nLin; i++)
        {
            tMean = 0;
            for (j = 0; j < nCol; j++)
                tMean += mMatrix(i,j);
            vMean.push_back(tMean/nCol);
        }
        vStd.clear();
        for (i = 0; i < nLin; i++)
        {
            tStd = 0;
            for (j = 0; j < nCol; j++)
                tStd += pow( (mMatrix(i,j) - vMean[i]),2);
            vStd.push_back(sqrt(tStd));
        }
        for (i = 0; i < nLin; i++)
        {
            for (j = 0; j < nLin; j++)
            {
                xyStd = 0;
                for(k = 0; k < nCol; k++)
                    xyStd += (mMatrix(i,k) - vMean[i])*(mMatrix(j,k) - vMean[j]);
                rMatrix(i,j) = xyStd/(vStd[i]*vStd[j]);
            }

        }
    }
    return rMatrix;
};
//---------------------------------------------------------------------
matrix<double> covariance(matrix<double> mMatrix, int nDimension)
{
    // nDimension == 0; correlation by lines => Sigma [cols x cols]
    // nDimension <> 0; correlation by cols => Sigma [lines x lines]
    int i, j, k, nLin, nCol;
    vector<double> vMean, vStd;
    double tMean, xyStd;
    nLin = mMatrix.RowNo();
    nCol = mMatrix.ColNo();
    matrix<double> rMatrix;

    if(nDimension == 0)
    {
        rMatrix.Null(nCol,nCol);

        vMean.clear();
        for (i = 0; i < nCol; i++)
        {
            tMean = 0.0;
            for (j = 0; j < nLin; j++)
                tMean += mMatrix(j,i);
            vMean.push_back(tMean/((float)nLin));
        }
        for (i = 0; i < nCol; i++)
        {
            for (j = i; j < nCol; j++)
            {
                xyStd = 0;
                for(k = 0; k < nLin; k++)
                    xyStd += (mMatrix(k,i) - vMean[i])*(mMatrix(k,j) - vMean[j]);
                rMatrix(i,j) = xyStd/((float)nLin - 1);
                rMatrix(j,i) = rMatrix(i,j);
            }

        }
    }
    else
    {
        rMatrix.Null(nLin,nLin);

        vMean.clear();
        for (i = 0; i < nLin; i++)
        {
            tMean = 0;
            for (j = 0; j < nCol; j++)
                tMean += mMatrix(i,j);
            vMean.push_back(tMean/((float)nCol));
        }
        for (i = 0; i < nLin; i++)
        {
            for (j = i; j < nLin; j++)
            {
                xyStd = 0.0;
                for(k = 0; k < nCol; k++)
                    xyStd += (mMatrix(i,k) - vMean[i])*(mMatrix(j,k) - vMean[j]);
                rMatrix(i,j) = xyStd/((float)nCol -1);
                rMatrix(j,i) = rMatrix(i,j);
            }
        }
    }
    return rMatrix;
};
//---------------------------------------------------------------------
matrix<double> StDNormalize(matrix<double> mMatrix, int nDimension, int afterIdx)
{
    // nDimension == 0; Normalize by lines
    // nDimension <> 0; normalize by cols
    // afterIdx = 0 normalize all data
    // afterIdx = 3 normalize from idx 3 to N-1 (init with 0)
    int i,j, nLin, nCol;
    vector<double> vMean, vStd;
    double tMean;
    nLin = mMatrix.RowNo();
    nCol = mMatrix.ColNo();
    matrix<double> mNormMatrix;

    mNormMatrix.SetSize(nLin,nCol);

    if(nDimension == 0)
    {
        vMean.clear();
        for (i = 0; i < nLin; i++)
        {
            tMean = 0;
            for (j = 0; j < nCol; j++)
                tMean += mMatrix(i,j);
            vMean.push_back(tMean/nCol);
        }
        vStd.clear();
        for (i = 0; i < nLin; i++)
        {
            tMean = 0;
            for (j = 0; j < nCol; j++)
                tMean += pow((mMatrix(i,j) - vMean[i]),2);
            vStd.push_back(sqrt(tMean/(nCol-1)));
        }
        for (i = 0; i < nLin; i++)
            for (j = 0; j < nCol; j++)
                if (i < afterIdx)
                    mNormMatrix(i,j) = mMatrix(i,j);
                else
                    mNormMatrix(i,j) = mMatrix(i,j)/vStd[i];
    }
    else
    {
        vMean.clear();
        for (i = 0; i < nCol; i++)
        {
            tMean = 0;
            for (j = 0; j < nLin; j++)
                tMean += mMatrix(j,i);
            vMean.push_back(tMean/nLin);
        }
        vStd.clear();
        for (i = 0; i < nCol; i++)
        {
            tMean = 0;
            for (j = 0; j < nLin; j++)
                tMean += pow((mMatrix(j,i) - vMean[i]),2);
            vStd.push_back(sqrt(tMean/(nLin-1)));
        }
        for (i = 0; i < nCol; i++)
            for (j = 0; j < nLin; j++)
                if (i < afterIdx)
                    mNormMatrix(j,i) = mMatrix(j,i);
                else
                    mNormMatrix(j,i) = mMatrix(j,i)/vStd[i];
    }
    return mNormMatrix;
};
//--------------------------------------------------------------------------
void swap(double &a, double &b)
{
    double tmp;
    tmp = a;
    a = b;
    b = tmp;
}
//--------------------------------------------------------------------------
void swap(int &a, int &b)
{
    int tmp;
    tmp = a;
    a = b;
    b = tmp;
}
//--------------------------------------------------------------------------
void swap(string &a, string &b)
{
    string tmp;
    tmp = a;
    a = b;
    b = tmp;
}
//--------------------------------------------------------------------------
int partition(vector<double> &vec, vector<int> &idx, int left, int right)
{
    int i, j;

    i = left;
    for (j = left + 1; j <= right; ++j)
    {
        if (vec[j] < vec[left])
        {
            ++i;
            swap(vec[i], vec[j]);
            swap(idx[i], idx[j]);
        }
    }
    swap(vec[left], vec[i]);
    swap(idx[left], idx[i]);

    return i;
}
//--------------------------------------------------------------------------
int partition(vector<string> &vec, vector<int> &idx, int left, int right)
{
    int i, j;

    i = left;
    for (j = left + 1; j <= right; ++j)
    {
        if (vec[j] < vec[left])
        {
            ++i;
            swap(vec[i], vec[j]);
            swap(idx[i], idx[j]);
        }
    }
    swap(vec[left], vec[i]);
    swap(idx[left], idx[i]);

    return i;
}
//--------------------------------------------------------------------------
void quickSort(vector<double> &vec, vector<int> &idx, int left, int right)
{
    int r;
    if (right > left)
    {
        r = partition(vec, idx, left, right);
        quickSort(vec, idx, left, r - 1);
        quickSort(vec, idx, r + 1, right);
    }
}
//--------------------------------------------------------------------------
void quickSort(vector<string> &vec, vector<int> &idx, int left, int right)
{
    int r;
    if (right > left)
    {
        r = partition(vec, idx, left, right);
        quickSort(vec, idx, left, r - 1);
        quickSort(vec, idx, r + 1, right);
    }
}
//--------------------------------------------------------------------------
void SortVector(vector<double> &vec, vector<int> &idx, int dec)
{
    // dec = 0 -> acend; dec = 1 -> descend
    int i;
    vector<double> tVec;
    vector<int> tIdx;

    if(idx.size() < 2)
    {
        idx.clear();
        for(i = 0; i < (int)vec.size(); i++)
            idx.push_back(i);
    }

    quickSort(vec, idx, 0, ((int)vec.size() - 1));
    if(dec == 1)
    {
        tVec = vec;
        tIdx = idx;
        idx.clear();
        vec.clear();
        for(i = ((int)tVec.size()-1); i >=0; i--)
        {
            vec.push_back(tVec[i]);
            idx.push_back(tIdx[i]);
        }
    }
};
//--------------------------------------------------------------------------
void SortVector(vector<string> &vec, vector<int> &idx, int dec)
{
    // dec = 0 -> acend; dec = 1 -> descend
    int i;
    vector<string> tVec;
    vector<int> tIdx;

    if(idx.size() < 2)
    {
        idx.clear();
        for(i = 0; i < (int)vec.size(); i++)
            idx.push_back(i);
    }

    quickSort(vec, idx, 0, ((int)vec.size() - 1));
    if(dec == 1)
    {
        tVec = vec;
        tIdx = idx;
        idx.clear();
        vec.clear();
        for(i = ((int)tVec.size()-1); i >=0; i--)
        {
            vec.push_back(tVec[i]);
            idx.push_back(tIdx[i]);
        }
    }
};
//--------------------------------------------------------------------------
void SortByVector(vector<string> &vec, vector<int> idx, int dec)
{
    // dec = 0 -> acend; dec = 1 -> descend
    int i;
    vector<string> tVec;
    tVec = vec;
    vec.clear();

    if(dec == 0)
        for(i = 0; i < (int)tVec.size(); i++)
            vec.push_back(tVec[idx[i]]);
    else
        for(i = ((int)tVec.size() -1); i >= 0; i--)
            vec.push_back(tVec[idx[i]]);
};
//--------------------------------------------------------------------------
vector<double> AcumFrequency(vector<double> mVector, int nBand)
{
    int i, n;
    double deltaN;
    vector<double> sVector;
    vector<int> idxVector;
    vector<double> fAcum;
    n = mVector.size();
    deltaN = n/(nBand+1);
    sVector = mVector;
    SortVector(sVector,idxVector,0);
    fAcum.clear();
    for(i = 0; i <= nBand; i++)
    {
        if(ceil(deltaN*i) >= n)
            fAcum.push_back(sVector[floor(deltaN*i)]);
        else
            fAcum.push_back(sVector[ceil(deltaN*i)]);
    }
    return fAcum;
};
//--------------------------------------------------------------------------
vector<double> OverSample(vector<double> mVector, int nSamples, int nSideLobs)
{
    vector<double> vOver, fSinc, vOut;
    int nPtSinc, i, k, xSize, idxIni, idxFim, xSizeOver, idxI, idxF;
    double xSinc;
    if (nSamples < 2)
        return mVector;
    else
        nSamples--;

    xSize = mVector.size();
    nPtSinc = 1 + 2*(nSideLobs*(nSamples + 1));
    xSizeOver = (xSize - 1)*(nSamples+1) + (nPtSinc - 1) + 1;
    for(i = 0; i < nPtSinc; i++)
    {
        xSinc = ((double)(i*2*M_PI*nSideLobs/(nPtSinc-1) - M_PI*nSideLobs));
        if (xSinc == 0)
            fSinc.push_back(1);
        else
            fSinc.push_back(sin(xSinc)/xSinc);
    }

    vOver.clear();
    for(i = 0; i < xSizeOver; i++)
        vOver.push_back(0);

    for(i = 0; i < xSize; i++)
    {
        idxI = i + i*nSamples;
        idxF = i + i*nSamples + nPtSinc;
        for(k = idxI; k < idxF; k++)
        {
            vOver[k] = vOver[k] + mVector[i]*fSinc[k-idxI];
        }
    }
    vOut.clear();
    idxIni = nSideLobs*(nSamples +1);
    idxFim = idxIni + (xSize - 1)*(nSamples+1) + 1;
    for(i = idxIni; i < idxFim; i++)
    {
        vOut.push_back(vOver[i]);
    };
    return vOut;
};
//--------------------------------------------------------------------------
vector<double> OverSample(vector<double> mVector, vector<double> &tVector, int nSamples, int nSideLobs)
{
    vector<double> vOver, fSinc, vOut;
    int nPtSinc, i, k, xSize, idxIni, idxFim, xSizeOver, idxI, idxF;
    double xSinc, tIni, tFim;
    if (nSamples < 2)
        return mVector;
    else
        nSamples--;

    xSize = mVector.size();
    nPtSinc = 1 + 2*(nSideLobs*(nSamples + 1));
    xSizeOver = (xSize - 1)*(nSamples+1) + (nPtSinc - 1) + 1;
    for(i = 0; i < nPtSinc; i++)
    {
        xSinc = ((double)(i*2*M_PI*nSideLobs/(nPtSinc-1) - M_PI*nSideLobs));
        if (xSinc == 0)
            fSinc.push_back(1);
        else
            fSinc.push_back(sin(xSinc)/xSinc);
    }


    vOver.clear();
    for(i = 0; i < xSizeOver; i++)
        vOver.push_back(0);

    for(i = 0; i < xSize; i++)
    {
        idxI = i + i*nSamples;
        idxF = i + i*nSamples + nPtSinc;
        for(k = idxI; k < idxF; k++)
        {
            vOver[k] = vOver[k] + mVector[i]*fSinc[k-idxI];
        }
    }
    vOut.clear();

    tIni = tVector[0];
    tFim = tVector[(tVector.size() -1)];
    tVector.clear();
    idxIni = nSideLobs*(nSamples +1);
    idxFim = idxIni + (xSize - 1)*(nSamples+1) + 1;

    for(i = idxIni; i < idxFim; i++)
    {
        tVector.push_back((tFim-tIni)*(i - idxIni)/(idxFim - idxIni-1) + tIni);
        vOut.push_back(vOver[i]);
    };
    return vOut;
};
//--------------------------------------------------------------------------
vector<double> FindPeaks(vector<double> vec, vector<int> &idx, int dec)
{
    // dec = 0 -> no sort; dec = 1 -> ascend; dec = 2 -> descend.
    int i, nSamples;
    vector<double> vPeaks;

    vPeaks.clear();
    idx.clear();

    nSamples = (int)vec.size();
    if (vec[0] > vec[1])
    {
        vPeaks.push_back(vec[0]);
        idx.push_back(0);
    }
    for(i = 1; i < (nSamples -1); i++)
        if ( (vec[i] > vec[i-1])  &&  (vec[i] > vec[i+1]) )
        {
            vPeaks.push_back(vec[i]);
            idx.push_back(i);
        }
    /*
    if (vec[(nSamples -1)] > vec[(nSamples -2)])
    {
        vPeaks.push_back(vec[(nSamples -1)]);
        idx.push_back((nSamples -1));
    }
    */
    if(dec == 0)
        return vPeaks;
    if(dec == 1)
    {
        SortVector(vPeaks,idx,0);
        return vPeaks;
    }
    if(dec == 2)
    {
        SortVector(vPeaks,idx,1);
        return vPeaks;
    }
    return vPeaks;
};
//--------------------------------------------------------------------------
matrix<double> GetRandonSample(matrix<double> mX, int iSig, double fError)
{
    matrix<double> Sample;
    double Zeta;
    int i, j, nSamples, nCols, nRows;
    vector<int> vIdx;

    if (iSig <= 90)
        Zeta = 1.65;
    if ((iSig > 90) || (iSig <= 95))
        Zeta = 1.96;
    if ((iSig > 95) || (iSig <= 97))
        Zeta = 2.17;
    if ((iSig > 97) || (iSig <= 98))
        Zeta = 2.33;
    if ((iSig > 99) || (iSig <= 99))
        Zeta = 2.58;
    if (iSig > 99)
        Zeta = 2.58;

    nCols = mX.ColNo();
    nRows = mX.RowNo();
    nSamples = pow((Zeta*0.5/fError),2);
    nSamples = ceil(nSamples/(1 + nSamples/nCols));

    for(i = 0; i < nCols; i++)
        vIdx.push_back(i);

    srand(unsigned( time(NULL) ));
    random_shuffle(vIdx.begin(), vIdx.end());
    Sample.SetSize(nRows,nSamples);
    for(i = 0; i < nSamples; i++)
        for(j = 0; j < nRows; j++)
            Sample(j,i) = mX(j,vIdx[i]);

    return Sample;

};
//--------------------------------------------------------------------------
matrix<double> csvread(string mFileName)
{
    FILE * pFile;
    int nLin, nCol, nColMax, nLinMax, c;
    string sValor;
    matrix<double> mX;
    mX.Null();
    char * pch;

    pFile = fopen (mFileName.c_str(),"r");
    nLin = 0;
    nCol = 0;
    nColMax = 1;
    nLinMax = 500;
    mX.SetSize(nLinMax,nColMax);
    sValor.clear();
    do {
        c = fgetc(pFile);
        sValor.append(1,c);
        if ((c == '\n') || (c == EOF))
        {
            pch = (char*)sValor.c_str();
            mX(nLin,nCol) = strtod(strtok(pch,","),NULL);
            sValor.clear();
            if (c != EOF)
            {
                nLin++;
                nCol = 0;
            };
            if (nLin == nLinMax)
            {
                nLinMax = nLinMax + 100;
                mX.SetSize(nLinMax,nColMax);
            };
        };
        if (c == ',')
        {
            pch = (char*)sValor.c_str();
            mX(nLin,nCol) = strtod(strtok(pch,","),NULL);
            sValor.clear();
            nCol++;
            if ((nCol == nColMax) && (nLin == 0))
            {
                nColMax++;
                mX.SetSize(nLinMax,nColMax);
            };
        }
    } while (c != EOF);
    mX.SetSize(nLin,nColMax);

    return mX;
};
//--------------------------------------------------------------------------
void csvwrite(matrix<double> mX, string mFileName)
{
    FILE * pFile;
    int i, j, nLin, nCol;
    string strTemp;
    char strValor[16];

    pFile = fopen (mFileName.c_str(),"w");

    nLin = mX.RowNo();
    nCol = mX.ColNo();
    for (i = 0; i < nLin; i++)
    {
        strTemp.clear();
        for(j = 0; j < nCol; j++)
        {
            sprintf(strValor,"%6.5e,",mX(i,j));
            strTemp.append(strValor);
        }
        strTemp.erase(strTemp.length()-1,1);
        fprintf (pFile,"%s\n",strTemp.c_str());
    }

    fclose (pFile);
};
//--------------------------------------------------------------------------
double InformacaoMutua(matrix<double> mX)
{
    int i, j, nLin;
    matrix<int> mA, mB;
    matrix<double> mX0, mX1, mX2, mX3;
    double X_2_3, X_2_15, sumA, sumB;
    double Imutual, sumRecImutual;

    nLin = mX.RowNo();
    if (nLin == 0)
        return 0;
    mA.SetSize(1,4);
    mB.SetSize(4,4);

    mX0 = NumSubElementos(mX,0, mA(0,0));
    mX1 = NumSubElementos(mX,1, mA(0,1));
    mX2 = NumSubElementos(mX,2, mA(0,2));
    mX3 = NumSubElementos(mX,3, mA(0,3));

    for(i = 0; i < 4; i++)
        NumSubElementos(mX0,i, mB(0,i));
    for(i = 0; i < 4; i++)
        NumSubElementos(mX1,i, mB(1,i));
    for(i = 0; i < 4; i++)
        NumSubElementos(mX2,i, mB(2,i));
    for(i = 0; i < 4; i++)
        NumSubElementos(mX3,i, mB(3,i));

    sumA = 0;
    sumB = 0;
    for(i = 0; i < 4; i++)
        sumA = pow(((double)mA(0,i) - (double)nLin/4),2);

    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            sumB = pow(((double)mB(j,i) - (double)nLin/16),2);

    X_2_3 = (16/9)*((double)1/nLin)*sumA;
    X_2_15 = (256/225)*((double)1/nLin)*sumB;
    if ((X_2_3 > 1.547) || (X_2_15 > 1.287))
    {
        sumRecImutual = InformacaoMutua(mX0) + InformacaoMutua(mX1) + InformacaoMutua(mX2) + InformacaoMutua(mX3);
        Imutual = nLin*log2(4) + sumRecImutual;
    }
    else
        Imutual = nLin*log2(nLin);

    return Imutual;
};
//--------------------------------------------------------------------------
//int NumSubElementos(matrix<double> mX, int iSec, matrix<double> &mY;)
matrix<double> NumSubElementos(matrix<double> mX, int iSec, int &NumElementos)
{
    int i, nLin, nCol, nHalfLin;
    vector<double> mX1, mX2, mOrdX1, mOrdX2;
    vector<int> mIdxX1, mIdxX2;
    matrix<double> mY;
    double medianaX1, medianaX2;
    mY.Null();
    nLin = mX.RowNo();
    nCol = mX.ColNo();
    NumElementos = 0;
    if (nLin == 0)
        return mY;

    mY.SetSize(nLin,nCol);
    mX1.clear();
    mX2.clear();
    for (i = 0; i < nLin; i++)
    {
       mX1.push_back (mX(i,0));
       mX2.push_back (mX(i,1));
    }

    switch (nLin)
    {
        case 1:
            medianaX1 = mX1[0];
            medianaX2 = mX2[0];
            break;
        case 2:
            medianaX1 = 0.5*(mX1[0] + mX1[1]);
            medianaX2 = 0.5*(mX2[0] + mX2[1]);
            break;
        default:
            mOrdX1 = mX1;
            mOrdX2 = mX2;
            SortVector(mOrdX1,mIdxX1,0);
            SortVector(mOrdX2,mIdxX2,0);
            if ((nLin%2) == 0)
            {
                nHalfLin = (int)(nLin/2);
                medianaX1 = 0.5*(mOrdX1[nHalfLin-1] + mOrdX1[nHalfLin]);
                medianaX2 = 0.5*(mOrdX2[nHalfLin-1] + mOrdX2[nHalfLin]);
            }
            else
            {
                nHalfLin = (int)((nLin + 1)/2);
                medianaX1 = 0.33333*(mOrdX1[nHalfLin-2] + mOrdX1[nHalfLin-1] + mOrdX1[nHalfLin]);
                medianaX2 = 0.33333*(mOrdX2[nHalfLin-2] + mOrdX2[nHalfLin-1] + mOrdX2[nHalfLin]);
            }
            break;
    }
    mY.SetSize(nLin,2);
    NumElementos = 0;
    switch (iSec)
    {
        case 0:
            for(i = 0; i < nLin; i++)
                if ( (mX(i,0) < medianaX1) && (mX(i,1) < medianaX2) )
                {
                    mY(NumElementos,0) = mX(i,0);
                    mY(NumElementos,1) = mX(i,1);
                    NumElementos++;
                };
            break;
        case 1:
            for(i = 0; i < nLin; i++)
                if ( (mX(i,0) < medianaX1) && (mX(i,1) >= medianaX2) )
                {
                    mY(NumElementos,0) = mX(i,0);
                    mY(NumElementos,1) = mX(i,1);
                    NumElementos++;
                };
            break;
        case 2:
            for(i = 0; i < nLin; i++)
                if ( (mX(i,0) >= medianaX1) && (mX(i,1) < medianaX2) )
                {
                    mY(NumElementos,0) = mX(i,0);
                    mY(NumElementos,1) = mX(i,1);
                    NumElementos++;
                };
            break;
        case 3:
            for(i = 0; i < nLin; i++)
                if ( (mX(i,0) >= medianaX1) && (mX(i,1) >= medianaX2) )
                {
                    mY(NumElementos,0) = mX(i,0);
                    mY(NumElementos,1) = mX(i,1);
                    NumElementos++;
                };
            break;
        default:
            break;
    }
    mY.SetSize(NumElementos,nCol);
    return mY;
};
//--------------------------------------------------------------------------
vector<double> spline(vector<double> X, vector<double> Y, vector<double> Xx)
{
    int k, i;
    int NkX, NkXx, Ord;
    double Ak, Bk, Ck, Dk, Yk;
    vector<double> Yy;
    matrix<double> Hk, Ek, Mk, Mtk;

    if ((X.size() != Y.size()) || (X.size() < 3))
            return Yy;

    NkX = (int)X.size();
    Ord = NkX - 1;
    NkXx = (int)Xx.size();

    Ek.Null();
    Hk.Null();
    Mtk.Null();
    Mk.Null();

    Ek.SetSize((Ord-1),1);
    Hk.SetSize((Ord-1),(Ord-1));
    Mtk.SetSize((Ord-1),1);
    Mk.SetSize(NkX,1);
    for(i = 0; i < (Ord-1); i++)
    {
        Ek(i,0) = (Y[i+2] - Y[i+1])/(X[i+2] - X[i+1]) - (Y[i+1] - Y[i])/(X[i+1] - X[i]);
        if (i == 0)
        {
            Hk(i,i) = (X[i+2] - X[i])/3;
            Hk(i,i+1) = (X[i+2] - X[i+1])/6;
            continue;
        }
        if (i == (Ord - 2))
        {
            Hk(i,i-1) = (X[i] - X[i-1])/6;
            Hk(i,i) = (X[i+2] - X[i])/3;
            break;
        }
        Hk(i,i-1) = (X[i] - X[i-1])/6;
        Hk(i,i) = (X[i+2] - X[i])/3;
        Hk(i,i+1) = (X[i+2] - X[i+1])/6;
    }
    Mtk = Hk.Solve(Ek);
    for(i = 0; i < NkX; i++)
    {
        if ((i == 0) || (i == (NkX-1)))
            Mk(i,0) = 0;
        else
            Mk(i,0) = Mtk(i-1,0);
    };

    k = 0;
    Yy.clear();
    for(i = 0; i < NkXx; i++)
    {
        if ( (Xx[i] >= X[k]) && (Xx[i] < X[k+1]) )
        {
            Ak = Y[k];
            Bk = (Y[k+1]  - Y[k])/(X[k+1]  - X[k]) - (Mk(k+1,0)/6 + Mk(k,0)/3)*(X[k+1]  - X[k]);
            Ck = Mk(k,0)/2;
            Dk = (Mk(k+1,0) - Mk(k,0))/(6*(X[k+1]  - X[k]));
            Yk = Ak + Bk*(Xx[i] - X[k]) + Ck*pow(Xx[i] - X[k],2) + Dk*pow(Xx[i] - X[k],3);
            Yy.push_back(Yk);
            continue;
        }
        if (k < (Ord-1))
        {
            k++;
            Ak = Y[k];
            Bk = (Y[k+1]  - Y[k])/(X[k+1]  - X[k]) - (Mk(k+1,0)/6 + Mk(k,0)/3)*(X[k+1]  - X[k]);
            Ck = Mk(k,0)/2;
            Dk = (Mk(k+1,0) - Mk(k,0))/(6*(X[k+1]  - X[k]));
            Yk = Ak + Bk*(Xx[i] - X[k]) + Ck*pow(Xx[i] - X[k],2) + Dk*pow(Xx[i] - X[k],3);
            Yy.push_back(Yk);
        }
    };
    return Yy;
};
//--------------------------------------------------------------------------
vector<string> read_directory(const string path, string sExt, bool dReadFileSize)
{
    vector<string> result;
    vector<int> idxFileListSize;
    vector<double> vFileSize;

    string temp, sFileName;
    unsigned int iLastPPos;
    struct stat file_status;
    dirent* de;
    DIR* dp;
    dp = opendir( path.empty() ? "." : path.c_str() );
    vFileSize.clear();
    if (dp)
    {
        while (true)
        {
            de = readdir( dp );
            if (de == NULL)
                break;
            if (sExt.size() == 0)
                result.push_back(string(de->d_name));
            else
            {
                if (de->d_type == DT_REG)
                {
                    sFileName = string(de->d_name);
                    iLastPPos = sFileName.find(sExt);

                    if ( (iLastPPos < sFileName.length()) && (iLastPPos == (sFileName.length() - sExt.length()))    )
                    {
                        result.push_back(string(de->d_name));
                        if (dReadFileSize)
                        {
                            temp.clear();
                            temp.append(path.c_str());
                            temp.append("/");
                            temp.append(de->d_name);
                            stat(temp.c_str(), &file_status);
                            vFileSize.push_back((double)file_status.st_size);
                        }
                    }
                }
                else
                    continue;

            }
        }
        closedir(dp);
        //sort( result.begin(), result.end() );
    }
    if (dReadFileSize)
    {
        SortVector(vFileSize,idxFileListSize,1); // Voltar para 1
        SortByVector(result,idxFileListSize,0);
    }
    else
         SortVector(result,idxFileListSize,0);

    return result;
};
//--------------------------------------------------------------------------
string RemoveExtension(string sNameFileExt)
{
    int pPos;
    pPos = sNameFileExt.find_last_of(".");
    return sNameFileExt.substr(0,pPos);
};
//--------------------------------------------------------------------------
vector<double> UnwrapFunc(vector<double> vX, bool bRemTrend)
{
    int i, k, nSamples;
    double Dp, Alpha, Beta, sumY, sumXY, sumX, sumX2;
    vector<double> vY;

    vY.clear();
    vY.push_back(vX[0]);

    nSamples = vX.size();
    for(i = 1; i < nSamples; i++)
        vY.push_back(0.0);

    k = 0;
    for(i = 1; i < nSamples; i++)
    {
        Dp = vX[i] - vX[i-1];
        if (Dp > M_PI)
            k--;
        if (Dp < -M_PI)
            k++;

        vY[i] = (vX[i] + 2*k*M_PI);

    }

    if (bRemTrend == true)
    {
        sumX = nSamples*(nSamples+1)/2;
        sumY = 0;
        sumXY = 0;
        sumX2 = 0;
        for(i = 0; i < nSamples; i++)
        {
            sumY += vY[i];
            sumXY += i*vY[i];
            sumX2 += i*i;
        };
        Alpha = (sumX2*sumY - sumXY*sumX)/((nSamples+1)*sumX2 - sumX*sumX);
        Beta = ((nSamples+1)*sumXY - sumX*sumY)/((nSamples+1)*sumX2 - sumX*sumX);;
        for(i = 0; i < nSamples; i++)
            vY[i] = vY[i] - Alpha - (double)i*Beta;
    }
    return vY;
};
//--------------------------------------------------------------------------
bool qr(matrix<double> A, matrix<double> &Q, matrix<double> &R)
{
    int i, j, k, n, nLin, nCol;
    double s, fak;//, mSum;
    matrix<double> I, P, T;
    vector<double> d;

    nLin = A.RowNo();
    nCol = A.ColNo();
    if (nLin != nCol)
        return false;

    I.Null(nLin,nCol);
    Q.Null(nLin,nCol);
    R.Null(nLin,nCol);
    for (i = 0; i < nLin; i++)
    {
        I(i,i) = 1.0;
        Q(i,i) = 1.0;
        d.push_back(0.0);
    }
    for(j = 0; j < nCol; j++)
    {
        s = 0.0;
        for (i = j; i < nLin; i++)
            s = s + pow(A(i,j),2);
        s = sqrt(s);

        if (A(j,j) > 0)
            d[j] = -s;
        else
            d[j] = s;

        fak = sqrt(s*(s+ abs(A(j,j))));

        A(j,j) = A(j,j) - d[j];
        for(k = j; k < nLin; k++)
            A(k,j) = A(k,j)/fak;

        for(i = j+1; i < nLin; i++)
        {
            s = 0.0;
            for(k = j; k < nLin; k++)
                s = s + A(k,j)*A(k,i);

            for(k = j; k < nLin; k++)
                A(k,i) = A(k,i) - A(k,j)*s;
        }

    }
    R.Null(nLin,nCol);
    for (i = 0; i < nLin; i++)
        for (j = i; j < nCol; j++)
        {
            if (i == j)
                R(i,j) = d[i];
            else
                R(i,j) = A(i,j);
        }
    //  conferir
    for(n = 0; n < nCol; n++)
    {
        P.Null(nLin,nCol);
        for(i = 0; i < nLin; i++)
        {
            for(j = 0; j < nCol; j++)
            {
                 if((i == j) && (i < n))
                    P(i,j) = 1;
                if((i == j) && (i >= n))
                    P(i,j) = 1 - A(i,n)*A(j,n);
                if((i != j) && (i >= n) && (j >= n))
                    P(i,j) = -A(i,n)*A(j,n);
            }

        }
        /*
        T.Null(nLin,nCol);
        for(i = 0; i < nLin; i++)
        {
            for(j = 0; j < nCol; j++)
            {
                mSum = 0.0;
                for(k = 0; k < nLin; k++)
                    mSum = mSum + Q(i,k)*P(k,j);
                T(i,j) = mSum;
            }
        }
        Q = T;
        */
        Q *= P;
    }
    return true;
};
//--------------------------------------------------------------------------
bool svd(matrix<double> A, matrix<double> &U, matrix<double> &S)
{
    int i, j, n, nLin, nCol, LoopMax, LoopCount;
    double dTol, dError, dNormE, dNormF;
    matrix<double> Q, R, Sn, St, E, F;

    nLin = A.RowNo();
    nCol = A.ColNo();
    if (nLin != nCol)
        return false;

    dTol=  2.2737e-13;//1024*DBL_EPSILON;
    dError = DBL_MAX;
    LoopMax = 100*nCol;
    LoopCount = 0;

    //I.Null(nLin,nCol);
    Q.Null(nLin,nCol);
    R.Null(nLin,nCol);
    U.Null(nLin,nCol);
    Sn.Null(nLin,nCol);
    St.Null(nLin,nCol);
    for (i = 0; i < nLin; i++)
    {
        U(i,i) = 1.0;
        for (j = 0; j < nCol; j++)
            Sn(i,j) = A(j,i);
    }
    while ( (dError > dTol) & (LoopCount < LoopMax) )
    {
        for (i = 0; i < nLin; i++)
            for (j = 0; j < nCol; j++)
                St(i,j) = Sn(j,i);
        qr(St,Q,Sn);
        U *= Q;

        for (i = 0; i < nLin; i++)
            for (j = 0; j < nCol; j++)
                St(i,j) = Sn(j,i);
        qr(St,Q,Sn);
        E.Null(nLin,nCol);
        F.Null(nLin,nCol);

        for (i = 0; i < nLin; i++)
        {
            F(i,i) = Sn(i,i);
            for (j = i+1; j < nCol; j++)
            {
                E(i,j) = Sn(i,j);
            }
        }
        dNormE = E.Norm();
        dNormF = F.Norm();
        if (dNormF == 0)
            dNormF = 1.0;
        dError = dNormE/dNormF;
        LoopCount++;
    }
    S.Null(nLin,nCol);
    for (n = 0; n < nLin; n++)
    {
        S(n,n) = abs(Sn(n,n));
        if (Sn(n,n) < 0)
            for(i = 0; i < nLin; i++)
                U(i,n) = -U(i,n);
    }
    return true;
};
//--------------------------------------------------------------------------
bool pca(matrix<double> A, matrix<double> &PCA, int nComp)
{
    int i, j, k, nDim, nSamples;
    //int L, C, nCol, nLin;
    double tPCA;
    matrix<double> corA, U, S;

    nDim = A.RowNo();
    nSamples = A.ColNo();
    corA = correlation(A,2);

    svd(corA,U,S);
    //nLin = corA.RowNo();
    //nCol = corA.ColNo();

    PCA.Null(nComp,nSamples);
    for(i = 0; i < nComp; i++)
        for(j = 0; j< nSamples; j++)
        {
            tPCA = 0;
            for(k = 0; k < nDim; k++)
                tPCA += U(i,k)*A(k,j);
            PCA(i,j) = tPCA;
        }
    return true;
};
//--------------------------------------------------------------------------
double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}
//--------------------------------------------------------------------------
