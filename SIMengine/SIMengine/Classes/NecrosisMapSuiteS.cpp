/*
 *  Created by Jan Poleszczuk
 *  Last modified 2020 by J Grajek
 */

#include "stdafx.h"
#include "NecrosisMapSuiteS.h"
#include <iostream>

NecrosisMap::NecrosisMap() {
    //setting default lattice sizes
    N1 = 100;
    N2 = 100;
    N3 = 1;
    numSpots = N1*N2*N3;
    afterInit = false;
};

NecrosisMap::~NecrosisMap() {
    if (afterInit) {
        delete[] changes;
        afterInit = false;
    }
};


void NecrosisMap::initialize(SIMparameters* parIn) {
    params = parIn;
    
    N1 = params->N1; N2 = params->N2; N3 = params->N3;
    numSpots = N1*N2*N3;
    
    N1D = N1 / 3; N2D = N2 / 3; N3D = N3/3;
    
    numSpotsD = N1D*N2D*N3D;
    
    
    if (afterInit) {
        delete[] changes;
    }
    afterInit = true;
    
    changes = DBG_NEW char[numSpots];
    memset(changes, 0, numSpots * sizeof(char));
    
    
    SpMat Snew(numSpotsD, numSpotsD);
    std::vector<T> coefficients;
    coefficients.reserve(7 * numSpotsD); //bo 7 przekatnych zapisujemy
    for (unsigned int i = 0; i < numSpotsD; ++i) {
        coefficients.push_back(T(i, i, params->DCnecro*6.));//0 diagonal
        if ((i + 1) % N1D)
            coefficients.push_back(T(i + 1, i, -params->DCnecro));//-1 diagonal
        if (i % N1D)
            coefficients.push_back(T(i - 1, i, -params->DCnecro));//1 diagonal
        if ((i % (N1D*N2D)) + N1D < N1D*N2D)
            coefficients.push_back(T(i + N1D, i, -params->DCnecro));//-N1 diagonal
        if (((int)(i % (N1D*N2D)) - (int)N1D) >= 0)
            coefficients.push_back(T(i - N1D, i, -params->DCnecro));//N1 diagonal
        if (i + N1D*N2D < numSpotsD)
            coefficients.push_back(T(i + N1D*N2D, i, -params->DCnecro));//-N1*N2 diagonal
        if (((int)i - (int)(N1D*N2D)) >= 0)
            coefficients.push_back(T(i - N1D*N2D, i, -params->DCnecro));//N1*N2 diagonal
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
                    b3D(k*N1D*N2D + j*N1D + i) = (double)numFS*params->DCnecro*params->physiologicalOxygen;
                    
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
                    bInterp3D(k*N1*N2 + j*N1 + i) = 1. / 3.*params->physiologicalOxygen;
                }
    
    state3D.resize(numSpots);
    state3D.fill(1.); //initial state
};

void NecrosisMap::addSink(unsigned int src, double glyc_rate) {
    changes[src] += 1;
}

void NecrosisMap::removeSink(unsigned int src, double glyc_rate) {
    changes[src] -= 1;
}

void NecrosisMap::updateMap() {
    
    char *iCC = DBG_NEW char[numSpotsD];
    memset(iCC, 0, numSpotsD*sizeof(char));
    int iC;
    int aux, cx, cy, cz;
    
    for (unsigned int i = 0; i < numSpots; ++i)
        if (changes[i] == 1) {
            aux = (int)(i % (N1*N2));
            cx = aux % (int)N1;
            cy = (int)((double)aux / (double)N1);
            cz = (int)((double)i / (double)(N1*N2));
            iC = (int)((cx / 3) + (cy / 3)*N1D + (cz / 3)*(N1D*N2D));
            iCC[iC]++;
        }
        else if (changes[i] == -1) {
            aux = (int)(i % (N1*N2));
            cx = aux % (int)N1;
            cy = (int)((double)aux / (double)N1);
            cz = (int)((double)i / (double)(N1*N2));
            iC = (cx / 3) + (cy / 3)*N1D + (cz / 3)*(N1D*N2D);
            iCC[iC]--;
        }
    
    int nnzU = 0, nnzD = 0;
    
    for (unsigned int i = 0; i < numSpotsD; i++)
        if (iCC[i]>0){
            nnzU++;
        }
        else if (iCC[i] < 0) {
            nnzD++;
        }
    if (nnzU || nnzD) {
        
        for (unsigned int i = 0; i < numSpotsD; ++i)
            if (iCC[i]!= 0) { //update
                RmS3D.coeffRef(i, i) += iCC[i]*params->TCnecro; // consumption where tumor cells are
            }
        
        solver3D.compute(RmS3D);
        if (solver3D.info() != Eigen::Success) {
            printf("NecrosisMap -> Compute step failed! \n");
            return;
        }
        Eigen::VectorXd x = solver3D.solve(b3D);
        if (solver3D.info() != Eigen::Success) {
            printf("Error: NecrosisMap -> Solving failed! \n");
            return;
        }
        //interpolating the state
        
        state3D = Interp3D*x + bInterp3D;
    }
    
    delete[]  iCC;
    //clearing changes vector
    memset(changes, 0, numSpots*sizeof(char));
}


mxArray* NecrosisMap::getState() {
    
    mwSize dims[3] = { (mwSize)N1, (mwSize)N2, (mwSize)N3 };
    mxArray *outMap = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *ptr = mxGetPr(outMap);
    
    memcpy(ptr, state3D.data(), numSpots*sizeof(double));
    
    
    return outMap;
}
