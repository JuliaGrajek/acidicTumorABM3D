/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "ChemotaxisMapSuiteS.h"
#include "NecrosisMapSuiteS.h"
#include "GlucoseMap.h"
#include "ProtonMap.h"
#include <iostream>
#include <stdio.h>
#include <tgmath.h> 

class Environment
{

private:

	struct point{//defining a point
		int x;
		int y;
		int z;
	};

	unsigned int *indcNeigh;//neighborhood
	bool *latticeOnlyTU;
	bool *necrosis;
    bool *hypoxia;
	bool *fibrosis;
    bool *latticeboundary;

	unsigned int numOccupiedSpots;
	unsigned int numNecroticSpots;
	unsigned int numFibroticSpots;

	float *IFNgMap;

	ChemotaxisMap ChtaxMap;
	NecrosisMap NecroMap;
    GlucoseMap GlucMap;
    ProtonMap ProtMap;
	

	unsigned int N1, N2, N3;
	unsigned int numSpots;

	CRandomMersenne* randGen; //pointer to random number generator
	SIMparameters* params; //pointer to parameters

	int numNeigh;
	unsigned int* neigh;
    unsigned int* neighLymph;
	float* chemoGrad;

	std::vector<point> IFNgMask, smoothMask;

public:
	unsigned int returnEmptyPlace(unsigned int);
	unsigned int returnEmptyPlaceImmune(unsigned int, float);
	
	unsigned int findTarget(unsigned int);

// 	float getAdjuValue(unsigned int pos){ return ChtaxMap.getValue(pos); };
	double getNecrosisValue(unsigned int pos){ return NecroMap.getValue(pos); };
    double getGlucoseValue(unsigned int pos){ return GlucMap.getValue(pos); };
    double getProtonValue(unsigned int pos){ return ProtMap.getValue(pos); };
    float getIFNgValue(unsigned int pos) { return IFNgMap[pos]; };

	Environment();
	~Environment();
	void initialize(SIMparameters*, CRandomMersenne*);
	void readState(const mxArray*);

	unsigned int getCenter();
	void newTUcell(unsigned int, double);
	void deleteTUcell(unsigned int, double);
    void newProtonSource(unsigned int pos, double gluc_uptake) {ProtMap.addSource(pos, gluc_uptake); };
    void clearProtonSources();
	void markNecrosis(unsigned int pos){ necrosis[pos] = true; };
    //void markHypoxia(unsigned int pos){ hypoxia[pos] = true; };
   //void unmarkHypoxia(unsigned int pos){ hypoxia[pos] = false; };

    void modulateHypoxia();
	void modulateIFNgMap(std::vector<unsigned int>::iterator, std::vector<unsigned int>::iterator, float);
	void modulateFibrosis(std::vector<unsigned int>::iterator, std::vector<unsigned int>::iterator);
	void decayIFNgMap();

	void updateNecroMap() { NecroMap.updateMap(); };
	void updateChemoMap() { ChtaxMap.updateMap(); };
    void updateGlucMap() { GlucMap.updateMap(); }
    void updateProtMap() { ProtMap.updateMap(); };
    

	unsigned int* generatePostions(int*);

	mxArray* getState();
};
