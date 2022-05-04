/*
*  Created by Jan Poleszczuk
*  Last modified May, 2022 by J Grajek
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "Environment.h"
#include <algorithm>
#include <iostream>

class TumorCells
{

public:
	struct TumorCell {//defining a tumor cell
		unsigned int place;
		unsigned char p;
		bool is_stem;
		float Antigen;
		unsigned char damage;
        double glyc_rate; //rate of glucose uptake, defaults to 1. If the cell enters a hypoxic region it adapts a higher glyc_rate (Pasteur effect);
        bool isAcidResistant;  // (e.g. expressing CA9?)
        bool PDL; 
	};

	void go_grow_die();
	
    double metabolize(TumorCell);

	void getStateForImmune(std::map<unsigned int, std::vector<TumorCells::TumorCell>::iterator>*);

	size_t NumTUcells();

	TumorCells(){};
	~TumorCells(){};

	void initialize(SIMparameters*, Environment*, CRandomMersenne*);
	void initializeFromState(const mxArray*, SIMparameters*, Environment*, CRandomMersenne*);
	void readState(mxArray*);

	mxArray* getState();

private:
	
    
	std::vector<TumorCell> cells; //vector containing all cells present in the system 
    std::vector<int> killedCells; // vector containing number of cells killed by T cells per round
	Environment* env; //pointer to the current environment
	CRandomMersenne* randGen; //pointer to random number generator
	SIMparameters* params; //pointer to parameters

};

