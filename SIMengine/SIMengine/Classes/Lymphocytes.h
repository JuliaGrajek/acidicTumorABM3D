/*
*  Created by Jan Poleszczuk
*  Last modified April, 2022 by J Grajek
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "Environment.h"
#include "TumorCells.h"
#include <iostream>

class Lymphocytes
{
private:
	struct Lymphocyte {//defining a lymphocyte
		unsigned int place;
		unsigned char p; //proliferating capacity
		unsigned char kcap; //killing capacity (how many tomor cells can be killed)
		unsigned char engaged; //flag telling whether the lymphocyte is attached to a tumor cell and can kill it
        bool quiescent; // flag telling whether the lymphocyte is in a place where pH is beneath the quiescence threshold
       
	};

	std::vector<Lymphocyte> cells; //vector containing all cells present in the system
	Environment* env; //pointer to the current environment
	CRandomMersenne* randGen; //pointer to random number generator
	SIMparameters* params; //pointer to parameters
	TumorCells* TUcells;

public:
	size_t NumLymphocytes();

	Lymphocytes(){};
	~Lymphocytes(){};
    
    std::vector<unsigned int> IFNgProducingCells; //vector containing information about Tcells that were close to Tumor cells and wanted to attack them (had found target). 
                                                  //In reality alla ctivated T cells produce IFNgamma, but we are only interested in the IFNgamma close to Tumor cells, so
                                                  //for computational piurposes we assume  only these T cells produce IFNgamma

	void initialize(SIMparameters*, Environment*, CRandomMersenne*, TumorCells*);
	void initializeFromState(const mxArray*, SIMparameters*, Environment*, CRandomMersenne*, TumorCells*);

	void influx();
	void influxInput(const mxArray*);

	void fibrosify();

	void action();
    
    void modulateIFNgMap();
    
	mxArray* getState();
};

