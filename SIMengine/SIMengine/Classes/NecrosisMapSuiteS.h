/*
 *  Created by Jan Poleszczuk
 *  Last modified 2020 by J Grajek
 */

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"

#include "suitesparse/cholmod.h"
//#include "suitesparse/cholmod_internal.h"


typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class NecrosisMap
{
    
private:
    
    unsigned int N1, N2, N3;
    unsigned int numSpots;
    
    unsigned int N1D, N2D, N3D;
    unsigned int numSpotsD;
    
    SIMparameters* params;

    
    SpMat RmS3D, Interp3D;
    Eigen::VectorXd b3D, bInterp3D, state3D;
    Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> solver3D;
    
    char* changes;
    bool afterInit;
    
public:
    
    NecrosisMap();
    ~NecrosisMap();
    
    void initialize(SIMparameters*);
    void updateMap();
    
    double getValue(unsigned int pos){
        return state3D(pos);
    };
    
    void addSink(unsigned int, double);
    void removeSink(unsigned int, double);
    
    mxArray* getState();
};
