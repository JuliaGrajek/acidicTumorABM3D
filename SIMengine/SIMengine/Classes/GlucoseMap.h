/*
 *  Created by Julia Grajek
 *  Last modified August, 2020 by Julia Grajek
 */

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"

#include "suitesparse/cholmod.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class GlucoseMap
{
    
private:
    
    unsigned int N1, N2, N3;
    unsigned int numSpots;
    
    unsigned int N1D, N2D, N3D;
    unsigned int numSpotsD;
    
    
    SIMparameters* params;
    
    cholmod_sparse *A, *Interp;
    cholmod_dense *state, *b, *bInterp;
    cholmod_factor *L;
    cholmod_common c;
    
    SpMat RmS3D, Interp3D;
    Eigen::VectorXd b3D, bInterp3D, state3D;
    Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> solver3D;
    
    
    double* changes_glyc_rate;
    bool afterInit;
    
public:
    
    GlucoseMap();
    ~GlucoseMap();
    
    void initialize(SIMparameters*);
    void updateMap();
    
    double getValue(unsigned int pos){
        return state3D(pos);
    };
    
    void addSink(unsigned int, double);
    void removeSink(unsigned int, double);
    
    mxArray* getState();
};
