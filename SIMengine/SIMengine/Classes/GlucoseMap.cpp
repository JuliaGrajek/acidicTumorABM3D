/*
 *  Created by Julia Grajek
 *  Last modified August, 2020 by Julia Grajek
 *  As of August 2020, we consider avascular tumor growth. Hence, we assume that blood vessels are located at the boundary of the domain.
 *  Glucose diffuses from there and is consumed by tumor cells
 */

#include "stdafx.h"
#include "GlucoseMap.h"
#include <iostream>


GlucoseMap::GlucoseMap() {
    //setting default lattice sizes
    N1 = 100;
    N2 = 100;
    N3 = 1;
    numSpots = N1*N2*N3;
    afterInit = false;
};

GlucoseMap::~GlucoseMap() {
    if (afterInit) {
        delete[] changes_glyc_rate;
        afterInit = false;
    }
};


void GlucoseMap::initialize(SIMparameters* parIn) {
    params = parIn;
    
    N1 = params->N1; N2 = params->N2; N3 = params->N3;
    numSpots = N1*N2*N3;
    
    N1D = N1 / 3; N2D = N2 / 3; N3D = N3 / 3;
    
    numSpotsD = N1D*N2D*N3D;
    
    if (afterInit) {
        delete[] changes_glyc_rate;
    }
    afterInit = true;
    
    changes_glyc_rate = DBG_NEW double[numSpots];
    memset(changes_glyc_rate, 0, numSpots*sizeof(double));
    
    SpMat Snew(numSpotsD, numSpotsD);
    std::vector<T> coefficients;
    coefficients.reserve(7 * numSpotsD); //bo 7 przekatnych zapisujemy
    for (unsigned int i = 0; i < numSpotsD; ++i) {
        coefficients.push_back(T(i, i, params->DCglucose*6.));//0 diagonal
        if ((i + 1) % N1D)
            coefficients.push_back(T(i + 1, i, -params->DCglucose));//-1 diagonal
        if (i % N1D)
            coefficients.push_back(T(i - 1, i, -params->DCglucose));//1 diagonal
        if ((i % (N1D*N2D)) + N1D < N1D*N2D)
            coefficients.push_back(T(i + N1D, i, -params->DCglucose));//-N1 diagonal
        if (((int)(i % (N1D*N2D)) - (int)N1D) >= 0)
            coefficients.push_back(T(i - N1D, i, -params->DCglucose));//N1 diagonal
        if (i + N1D*N2D < numSpotsD)
            coefficients.push_back(T(i + N1D*N2D, i, -params->DCglucose));//-N1*N2 diagonal
        if (((int)i - (int)(N1D*N2D)) >= 0)
            coefficients.push_back(T(i - N1D*N2D, i, -params->DCglucose));//N1*N2 diagonal
    }
    
    
    Snew.setFromTriplets(coefficients.begin(), coefficients.end());
    Snew.makeCompressed();
    RmS3D = Snew;
    
    b3D.resize(numSpotsD);
    b3D.fill(0.);
    //vector b
    for (unsigned int i = 0; i < N1D; ++i)
        for (unsigned int j = 0; j < N2D; ++j)
            for (unsigned int k = 0; k < N3D; ++k)
                if (i == 0 || i == N1D - 1 || j == 0 || j == N2D - 1 || k == 0 || k == N3D - 1) {
                    char numFS = (char)(i == 0) + (char)(i == N1D - 1) + (char)(j == 0) + (char)(j == N2D - 1) + (char)(k == 0) + (char)(k == N3D - 1);
                    b3D(k*N1D*N2D + j*N1D + i) = (double)numFS*params->DCglucose*params->physiologicalGlucose;
                }
    
    
// 		initializing Interpolation matrix
    coefficients.clear();
    SpMat InterpTmp(numSpots, numSpotsD);
    int aux, cx, cy, cz, whX, whY, whZ, rPx, rPy, rPz, whX2, whY2, whZ2;
    for (unsigned int i = 0; i < numSpots; ++i) {
        //casting to super node
        aux = (int)(i % (N1*N2));
        cx = aux % (int)N1;
        cy = (int)((double)aux / (double)N1);
        cz = (int)((double)i / (double)(N1*N2));
        
        whX = cx / 3;
        whY = cy / 3;
        whZ = cz / 3;
        
        rPx = (cx % 3) - 1;
        rPy = (cy % 3) - 1;
        rPz = (cz % 3) - 1;
        
        if (rPx != 0 || rPy != 0 || rPz != 0) {
            whX2 = whX + rPx;
            whY2 = whY + rPy;
            whZ2 = whZ + rPz;
            
            coefficients.push_back(T(i, whX + whY*N1D + whZ*N1D*N2D, 2. / 3.));//closer super node that is in the same sparse grid
            
            if (whX2 >= 0 && whX2 < (int)N1D && whY2 >= 0 && whY2 < (int)N2D && whZ2 >= 0 && whZ2 < (int)N3D)
                coefficients.push_back(T(i, whX2 + whY2*N1D + whZ2*N1D*N2D, 1. / 3.));//other node
        }
        else {
            coefficients.push_back(T(i, whX + whY*N1D + whZ*N1D*N2D, 1.));//center super node
        }
    }
    
    InterpTmp.setFromTriplets(coefficients.begin(), coefficients.end());
    InterpTmp.makeCompressed();
    Interp3D = InterpTmp;
    
    bInterp3D.resize(numSpots);
    bInterp3D.fill(0.);
    
    for (unsigned int i = 0; i < N1; ++i)
        for (unsigned int j = 0; j < N2; ++j)
            for (unsigned int k = 0; k < N3; ++k)
                if (i == 0 || i == N1 - 1 || j == 0 || j == N2 - 1 || k == 0 || k == N3 - 1) {
                    bInterp3D(k*N1*N2 + j*N1 + i) = 1. / 3.*params->physiologicalGlucose;
                }
    
    state3D.resize(numSpots);
    state3D.fill(1.); //initial state
    
    
};

void GlucoseMap::addSink(unsigned int src, double gr) {
    changes_glyc_rate[src] +=1;
}

void GlucoseMap::removeSink(unsigned int src, double gr) {
    changes_glyc_rate[src] -= 1;
}

void GlucoseMap::updateMap() {
    
    double *iCC = DBG_NEW double[numSpotsD];
    memset(iCC, 0, numSpotsD*sizeof(double));
    
    int iC;
    int aux, cx, cy, cz;
    
    for (unsigned int i = 0; i < numSpots; ++i)
        if (changes_glyc_rate[i] ==1) {
            aux = (int)(i % (N1*N2));
            cx = aux % (int)N1;
            cy = (int)((double)aux / (double)N1);
            cz = (int)((double)i / (double)(N1*N2));
            iC = (int)((cx / 3) + (cy / 3)*N1D + (cz / 3)*(N1D*N2D));
            iCC[iC]+=1;
        }
        else if (changes_glyc_rate[i] ==-1) {
            aux = (int)(i % (N1*N2));
            cx = aux % (int)N1;
            cy = (int)((double)aux / (double)N1);
            cz = (int)((double)i / (double)(N1*N2));
            iC = (cx / 3) + (cy / 3)*N1D + (cz / 3)*(N1D*N2D);
            iCC[iC]-=1;
        }
    
    int nnzU = 0, nnzD = 0; // nnzUL = 0, nnzDL = 0;
    
    for (unsigned int i = 0; i < numSpotsD; i++)
        if (iCC[i]>0){
            nnzU++;
        }
        else if (iCC[i] < 0) {
            nnzD++;
        }
    
    
    
    if (nnzU || nnzD ){
        
        for (unsigned int i = 0; i < numSpotsD; ++i)
            if (iCC[i]!= 0) { //update
                RmS3D.coeffRef(i, i) += iCC[i]*params->TCglucose; // consumption where tumor cells are
            }
        
        
        solver3D.compute(RmS3D);
        if (solver3D.info() != Eigen::Success) {
            printf("GlucoseMap -> Compute step failed! \n");
            return;
        }
        Eigen::VectorXd x = solver3D.solve(b3D);
        if (solver3D.info() != Eigen::Success) {
            printf("Error: GlucoseMap -> Solving failed! \n");
            std::cout << x;
            return;
        }
        //interpolating the state
        state3D = Interp3D*x + bInterp3D;
        
    }
    
    
    delete[] iCC;
    memset(changes_glyc_rate, 0, numSpots*sizeof(double));
    
}


mxArray* GlucoseMap::getState() {
    
    mwSize dims[3] = { (mwSize)N1, (mwSize)N2, (mwSize)N3 };
    mxArray *outMap = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *ptr = mxGetPr(outMap);
    
    memcpy(ptr, state3D.data(), numSpots*sizeof(double));
    
    
    return outMap;
}