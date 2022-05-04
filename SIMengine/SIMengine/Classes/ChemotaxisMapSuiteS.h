/*
 *  Created by Jan Poleszczuk
 *  Last modified April 2022, by Julia Grajek
 */

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "suitesparse/cholmod.h"


typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class ChemotaxisMap
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
    
    ChemotaxisMap();
    ~ChemotaxisMap();
    
    void initialize(SIMparameters*);
    void updateMap();
    
    float getValue(unsigned int pos){
        
        return (float)state3D(pos);
        
    };
    
    void addSource(unsigned int);
    void removeSource(unsigned int);
    
    mxArray* getState();
    
};
