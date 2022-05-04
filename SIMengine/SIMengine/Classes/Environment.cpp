/*
 *  Created by Jan Poleszczuk
 *  Last modified April, 2022 by J Grajek
 */

#include "stdafx.h"
#include "Environment.h"

Environment::Environment() {
    latticeboundary = NULL;
    indcNeigh = NULL;
    neigh = NULL;
    neighLymph = NULL;
    necrosis = NULL;
    hypoxia = NULL;
    fibrosis = NULL;
    IFNgMap = NULL;
    latticeOnlyTU = NULL;
    chemoGrad = NULL;
    //setting default lattice sizes
    N1 = 10;
    N2 = 10;
    N3 = 10;
    numSpots = 1000;
    numNeigh = 26;
    numOccupiedSpots = 0;
    numNecroticSpots = 0;
    numFibroticSpots = 0;
};

Environment::~Environment() {
    if (latticeboundary != NULL)
        delete[] latticeboundary;
    if (indcNeigh != NULL)
        delete[] indcNeigh;
    if (neigh != NULL)
        delete[] neigh;
    if (neighLymph != NULL)
        delete[] neighLymph;
    if (necrosis != NULL)
        delete[] necrosis;
    if (hypoxia != NULL)
        delete[] hypoxia;
    if (fibrosis != NULL)
        delete[] fibrosis;
    if (IFNgMap != NULL)
        delete[] IFNgMap;
    if (latticeOnlyTU != NULL)
        delete[] latticeOnlyTU;
    if (chemoGrad != NULL)
        delete[] chemoGrad;
}

unsigned int Environment::getCenter() {
    return N1 / 2 + (N2 / 2)*N1 + (N3 / 2)*N1*N2; // integer division!!!
}


void Environment::newTUcell(unsigned int pos, double glyc_rate) {
    latticeOnlyTU[pos] = true;
    necrosis[pos] = false;
    numOccupiedSpots++;
    ChtaxMap.addSource(pos);
    NecroMap.addSink(pos, glyc_rate);
    GlucMap.addSink(pos, glyc_rate);
}

void Environment::deleteTUcell(unsigned int pos, double glyc_rate) {
    latticeOnlyTU[pos] = false;
    numOccupiedSpots--;
    ChtaxMap.removeSource(pos);
    NecroMap.removeSink(pos, glyc_rate);
    GlucMap.removeSink(pos, glyc_rate);
}

void Environment::clearProtonSources() {
    for (unsigned int i = 0; i < numSpots; ++i){
        ProtMap.removeSource(i);
    }
}

unsigned int* Environment::generatePostions(int* N) {
    
    unsigned int* positions = NULL;
    
    std::vector<unsigned int> availableSpots;
    
    availableSpots.reserve(numSpots - numOccupiedSpots - numNecroticSpots);
    
    
    int r = (int) (std::min({N1, N2, N3})/2);
    
    for (unsigned int i = 0; i < numSpots; ++i){
        int aux = (int)((i) % (N1*N2));
        int cx = (int)(N1/2)-(aux % (int)N1);
        int cy = (int)(N2/2)- (int)((double)aux / (double)N1);
        int cz = (int)(N3/2)-(int)((double)(i) / (double)(N1*N2));
        
        if ( (pow(cx,2)+ pow(cy,2)+pow(cz,2)<pow(r,2)) && !latticeOnlyTU[i] &&!latticeboundary[i] && !necrosis[i] && (!fibrosis[i] || randGen->Random() < params->stromaPerm))
            availableSpots.push_back(i);}
    
    if (availableSpots.size() < *N) {
        *N = (int)availableSpots.size();
        if (*N > 0) {
            
            positions = DBG_NEW unsigned int[*N];
            for (int i = 0; i < *N; ++i)
                positions[i] = availableSpots.at(i);
            
        }
    }
    else {
        
        positions = DBG_NEW unsigned int[*N];
        int gP;
        for (int i = 0; i < *N; ++i) {
            gP = randGen->Rand((int)availableSpots.size());
            positions[i] = availableSpots.at(gP);
            availableSpots.erase(availableSpots.begin() + gP);
            
        }
    }
    availableSpots.clear();
    return positions;
}

void Environment::decayIFNgMap() {
    for (unsigned int i = 0; i < numSpots; ++i)
        if (IFNgMap[i] > 0)
            IFNgMap[i] *= params->IFNgDecay;
}

void Environment::modulateFibrosis(std::vector<unsigned int>::iterator bg, std::vector<unsigned int>::iterator ed) {
    for (std::vector<unsigned int>::iterator it = bg; it != ed; it++) {//iterate through seeds
        fibrosis[(*it)] = true;
        
        //getting coordinates
        int aux = (int)((*it) % (N1*N2));
        int cx = aux % (int)N1;
        int cy = (int)((double)aux / (double)N1);
        int cz = (int)((double)(*it) / (double)(N1*N2));
        
        //going through neighborhood
        unsigned int pos;
        for (int i = 0; i < smoothMask.size(); ++i) {
            int cxN = cx + smoothMask.at(i).x;
            int cyN = cy + smoothMask.at(i).y;
            int czN = cz + smoothMask.at(i).z;
            if (cxN >= 0 && cxN < (int)N1 && cyN >= 0 && cyN < (int)N2 && czN >= 0 && czN < (int)N3 && randGen->Random() < params->fibrFrac) {
                pos = (unsigned int)czN*N1*N2 + (unsigned int)cyN*N1 + (unsigned int)cxN;
                if (params->seedUnderneath) {
                    fibrosis[pos] = true;
                } else if(!latticeOnlyTU[pos] && !latticeboundary[pos] ) {//if the lattice spot is empty
                    fibrosis[pos] = true;
                }
            }
        }
    }
}

void Environment::modulateIFNgMap(std::vector<unsigned int>::iterator bg, std::vector<unsigned int>::iterator ed, float addIFNg) {
    
    for (std::vector<unsigned int>::iterator it = bg; it != ed; it++) {//iterate through seeds
        //getting coordinates
        int aux = (int)((*it) % (N1*N2));
        int cx = aux % (int)N1;
        int cy = (int)((double)aux / (double)N1);
        int cz = (int)((double)(*it) / (double)(N1*N2));
        
        //going through neighborhood
        unsigned int pos;
        for (int i = 0; i < IFNgMask.size(); ++i) {
            int cxN = cx + IFNgMask.at(i).x;
            int cyN = cy + IFNgMask.at(i).y;
            int czN = cz + IFNgMask.at(i).z;
            if (cxN >= 0 && cxN < (int)N1 && cyN >= 0 && cyN < (int)N2 && czN >= 0 && czN < (int)N3) {
                pos = (unsigned int)czN*N1*N2 + (unsigned int)cyN*N1 + (unsigned int)cxN;
                IFNgMap[pos] += addIFNg;
            }
        }
    }
}

void Environment::modulateHypoxia(){
    numSpots = params->N1*params->N2*params->N3;
    for (unsigned int pos = 0; pos < numSpots; ++pos) {//iterate through seeds
        if (NecroMap.getValue(pos)< params->hypThresh && necrosis[pos]==false){
            hypoxia[pos] = true;
        }
        else {
            hypoxia[pos] = false;
        }
    }
}


unsigned int Environment::returnEmptyPlaceImmune(unsigned int indx, float CErwalk) {
    int nF = 0;
    float maxChemo = 0.0f;
    
    for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
        if (!latticeboundary[indx - indcNeigh[j]] && ( !fibrosis[indx - indcNeigh[j]]  || randGen->Random() < params->stromaPerm)) {
            neigh[nF] = indx - indcNeigh[j];
            chemoGrad[nF] = ChtaxMap.getValue(neigh[nF]);
            if (chemoGrad[nF] > maxChemo)
                maxChemo = chemoGrad[nF];
            nF++;
        }
    }
    for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
        if (!latticeboundary[indx + indcNeigh[j]] && ( !fibrosis[indx + indcNeigh[j]]  || randGen->Random() < params->stromaPerm)) {
            neigh[nF] = indx + indcNeigh[j];
            chemoGrad[nF] = ChtaxMap.getValue(neigh[nF]);
            if (chemoGrad[nF] > maxChemo)
                maxChemo = chemoGrad[nF];
            nF++;
        }
    }
    if (nF) {//selecting free spot taking into account chemotaxis
        
        if (maxChemo > 0){
            float minChemo = (1.0f - CErwalk)*chemoGrad[0] / maxChemo + CErwalk*(float)randGen->Random();
            int which = 0;
            for (int j = 1; j < nF; j++) {
                float propMinChemo = (1.0f - CErwalk)*chemoGrad[j] / maxChemo + CErwalk*(float)randGen->Random();
                if (propMinChemo > minChemo) {
                    which = j;
                    minChemo = propMinChemo;
                }
            }
            
            return neigh[which];
        }
        else {
            return neigh[randGen->Rand(nF)];
        }
    }
    else {//no free spot
        return 0;
    }
}

unsigned int Environment::returnEmptyPlace(unsigned int indx) {
    int nF = 0;
    for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
        if (!latticeOnlyTU[indx - indcNeigh[j]] && !latticeboundary[indx-indcNeigh[j]]  && (!fibrosis[indx - indcNeigh[j]] || randGen->Random() < params->stromaPerm) ) {
            neigh[nF] = indx - indcNeigh[j];
            nF++;
        }
    }
    for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
        if (!latticeOnlyTU[indx + indcNeigh[j]]  && !latticeboundary[indx+indcNeigh[j]]  && (!fibrosis[indx + indcNeigh[j]] || randGen->Random() < params->stromaPerm)) {
            neigh[nF] = indx + indcNeigh[j];
            nF++;
        }
    }
    if (nF) {//selecting free spot at random
        return neigh[randGen->Rand(nF)];
    }
    else {//no free spot
        return 0;
    }
}

unsigned int Environment::findTarget(unsigned int indx) {
    int nF = 0;
    for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
        if (latticeOnlyTU[indx - indcNeigh[j]]) {
            neighLymph[nF] = indx - indcNeigh[j];
            nF++;
        }
    }
    for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
        if (latticeOnlyTU[indx + indcNeigh[j]]) {
            neighLymph[nF] = indx + indcNeigh[j];
            nF++;
        }
    }
    if (latticeOnlyTU[indx]){
        neighLymph[nF] = indx;
        nF++;
    }
    
    if (nF) {//selecting free spot at random
        
        return neighLymph[randGen->Rand(nF)];
    }
    else {//no free spot
        return 0;
    }
}

void Environment::initialize(SIMparameters* parIn, CRandomMersenne *rG) {
    //setting enviroment dimension
    params = parIn;
    N1 = params->N1; N2 = params->N2; N3 = params->N3;
    numSpots = N1*N2*N3;
    
    numOccupiedSpots = 0;
    numNecroticSpots = 0;
    numFibroticSpots = 0;
    
    randGen = rG;
    
    if (necrosis != NULL)
        delete[] necrosis;
    if (hypoxia != NULL)
        delete[] hypoxia;
    if (fibrosis != NULL)
        delete[] fibrosis;
    if (IFNgMap != NULL)
        delete[] IFNgMap;
    if (latticeOnlyTU != NULL)
        delete[] latticeOnlyTU;
    if (latticeboundary != NULL)
        delete[] latticeboundary;
    
    IFNgMask.clear();
    smoothMask.clear();
    
    necrosis = DBG_NEW bool[numSpots];
    memset(necrosis, 0, numSpots * sizeof(bool));
    hypoxia = DBG_NEW bool[numSpots];
    memset(hypoxia, 0, numSpots * sizeof(bool));
    fibrosis = DBG_NEW bool[numSpots];
    memset(fibrosis, 0, numSpots * sizeof(bool));
    latticeOnlyTU = DBG_NEW bool[numSpots];
    memset(latticeOnlyTU, 0, numSpots * sizeof(bool));
    latticeboundary = DBG_NEW bool[numSpots];
    memset(latticeboundary, 0, numSpots * sizeof(bool));
    IFNgMap = DBG_NEW float[numSpots];
    memset(IFNgMap, 0, numSpots * sizeof(float));
    
    ChtaxMap.initialize(params);
    NecroMap.initialize(params);
    GlucMap.initialize(params);
    ProtMap.initialize(params);
    
    //setting lattice boundary
    for (unsigned int i = 0; i < N1; ++i)
        for (unsigned int j = 0; j < N2; ++j)
            for (unsigned int k = 0; k < N3; ++k)
                if (i == 0 || i == N1 - 1 || j == 0 || j == N2 - 1 || k == 0 || k == N3 - 1) {
                    latticeboundary[k*N1*N2 + j*N1 + i] = true;
                    numOccupiedSpots++;
                }
    
    //defining the neighborhood
    if (indcNeigh != NULL)
        delete[] indcNeigh;
    indcNeigh = DBG_NEW unsigned int[13];
    int indx, aux = 0;
    for (int i = -1; i < 2; ++i)
        for (int j = -1; j < 2; ++j)
            for (int k = -1; k < 2; ++k) {
                indx = k*N1*N2+j*N1+i;
                if (indx > 0) {
                    indcNeigh[aux] = (unsigned int)indx;
                    aux++;
                }
                
            }
    
    if (neigh != NULL)
        delete[] neigh;
    neigh = DBG_NEW unsigned int[26];
    numNeigh = 13;
    
    if (neighLymph != NULL)
        delete[] neighLymph;
    neighLymph = DBG_NEW unsigned int[27];
    
    
    if (chemoGrad != NULL)
        delete[] chemoGrad;
    chemoGrad = DBG_NEW float[26];
    //defining IFNg mask
    int r2 = params->IFNgRange*params->IFNgRange;
    for (int i = -params->IFNgRange; i <= params->IFNgRange; ++i)
        for (int j = -params->IFNgRange; j <= params->IFNgRange; ++j)
            for (int k = -params->IFNgRange; k <= params->IFNgRange; ++k)
                if (i*i + j*j + k*k <= r2) {//within the sphere
                    point point = { i, j, k };
                    IFNgMask.push_back(point);
                }
    //defining smoothMask
    r2 = params->smoothRadius*params->smoothRadius;
    for (int i = -params->smoothRadius; i <= params->smoothRadius; ++i)
        for (int j = -params->smoothRadius; j <= params->smoothRadius; ++j)
            for (int k = -params->smoothRadius; k <= params->smoothRadius; ++k)
                if (i*i + j*j + k*k <= r2) {//within the sphere
                    point point = { i, j, k };
                    smoothMask.push_back(point);
                }
}

void Environment::readState(const mxArray* sysTempl) {
    mxArray *grid = mxGetField(sysTempl, 0, "grid");
    float* IFNgMapS = (float*)mxGetPr(mxGetField(grid, 0, "IFNgMap"));
    memcpy(IFNgMap, IFNgMapS, numSpots * sizeof(float));
    bool* Lf = (bool*)mxGetPr(mxGetField(grid, 0, "Lf"));
    memcpy(fibrosis, Lf, numSpots * sizeof(bool));
    
}

mxArray* Environment::getState() {
    
    const char *field_names[] = { "Lb", "Ln", "Lh", "Lf", "IFNgMap", "ChtaxMap","NecroMap", "GlucMap", "ProtMap"}; // L - matrix with occupied spaces, Ln - matrix with places with necrosis, Lf - matrix with fibrosis, Lh - matrix with hypoxia
    mwSize dimsStruct[2] = { 1, 1 };
    mxArray* envState = mxCreateStructArray(1, dimsStruct, 9, field_names);
    
    //2 - copying the L matrix
    mwSize dims[3] = { (mwSize)N1, (mwSize)N2, (mwSize)N3 };
    mxSetField(envState, 0, "Lb", mxCreateLogicalArray(3, dims));
    bool *Lb = (bool*)mxGetPr(mxGetField(envState, 0, "Lb"));
    memcpy(Lb, latticeboundary, numSpots * sizeof(bool));
    
    mxSetField(envState, 0, "Ln", mxCreateLogicalArray(3, dims));
    bool *Ln = (bool*)mxGetPr(mxGetField(envState, 0, "Ln"));
    memcpy(Ln, necrosis, numSpots * sizeof(bool));
    
    mxSetField(envState, 0, "Lh", mxCreateLogicalArray(3, dims));
    bool *Lh = (bool*)mxGetPr(mxGetField(envState, 0, "Lh"));
    memcpy(Lh, hypoxia, numSpots * sizeof(bool));
    
    mxSetField(envState, 0, "Lf", mxCreateLogicalArray(3, dims));
    bool *Lf = (bool*)mxGetPr(mxGetField(envState, 0, "Lf"));
    memcpy(Lf, fibrosis, numSpots * sizeof(bool));
    
    mxSetField(envState, 0, "IFNgMap", mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL));
    float *Imap = (float*)mxGetPr(mxGetField(envState, 0, "IFNgMap"));
    memcpy(Imap, IFNgMap, numSpots * sizeof(float));
    
    mxSetField(envState, 0, "ChtaxMap", ChtaxMap.getState());
    mxSetField(envState, 0, "NecroMap", NecroMap.getState());
    mxSetField(envState, 0, "GlucMap", GlucMap.getState());
    mxSetField(envState, 0, "ProtMap", ProtMap.getState());
    return envState;
    
}
