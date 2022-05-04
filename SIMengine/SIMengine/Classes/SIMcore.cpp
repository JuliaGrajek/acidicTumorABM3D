/*
*  Created by Jan Poleszczuk
*  Last modified September, 2017 by Jan Poleszczuk
*/


#include "SIMcore.h"

SIMcore::SIMcore(){};
SIMcore::~SIMcore(){};

void SIMcore::initialize(const mxArray* sysTempl, const mxArray* cnst) {

	//initializing parameters
	params.initialize(sysTempl, cnst);
	//params.listSIMParameters();
    
	//intitializing random number generator
	if (mxIsNaN(params.initSeed)) {//no seed specified
        //mexPrintf("No initial seed specified!\n");
		randGen.RandomInit(time(NULL));
	}
	else {
        //mexPrintf("Initial seed specified!\n");
		randGen.RandomInit((int)params.initSeed);
		randGen.LastInterval = 0;
	}

	//initializing environment (important to be before agents initialization)
	env.initialize(&params, &randGen);

	//initializing cells
	TUcells.initialize(&params, &env, &randGen);
	lymphocytes.initialize(&params, &env, &randGen, &TUcells);
}

void SIMcore::initializeMexFree() {

	params.listSIMParameters();

	//intitializing random number generator
	randGen.RandomInit((int)params.initSeed);
	randGen.LastInterval = 0;

	//initializing environment
	env.initialize(&params, &randGen);

	//initializing cells
	TUcells.initialize(&params, &env, &randGen);
	lymphocytes.initialize(&params, &env, &randGen, &TUcells);
}

void SIMcore::initializeFromState(const mxArray* sysTempl, const mxArray* cnst) {

	//initializing parameters
	params.initialize(sysTempl, cnst);
	//params.listSIMParameters();

	//intitializing random number generator
     if (mxIsNaN(params.initSeed)) {//no seed specified
        //mexPrintf("No initial seed specified!\n");
        randGen.RandomInit(time(NULL));
    }
    else {
        //mexPrintf("Initial seed specified!\n");
        randGen.RandomInit((int)params.initSeed);
        randGen.LastInterval = 0;
    }

	//initializing environment (important to be before agents initialization)
	env.initialize(&params, &randGen);
	env.readState(sysTempl);

	//initializing cells
	TUcells.initializeFromState(sysTempl, &params, &env, &randGen);
    double * exp_id;
    double * lym_id;
    exp_id = mxGetPr(mxGetField(sysTempl,0, "experiment_id"));
    lym_id = mxGetPr(mxGetField(sysTempl,0, "lym_id"));
//     std::cout << "lym_id: " << lym_id << '\n';
    if (*exp_id==1 || *lym_id==1){
        lymphocytes.initializeFromState(sysTempl, &params, &env, &randGen, &TUcells);
        
    }
    else{
        lymphocytes.initialize(&params, &env, &randGen, &TUcells); 
    }
	//updating maps
	env.updateChemoMap();
	env.updateNecroMap();
    env.updateGlucMap();
    env.updateProtMap();
    
}

void SIMcore::readState(const mxArray *state) {
	mxArray *TUc = mxGetField(state, 0, "TUcells");
	TUcells.readState(TUc);
}



size_t SIMcore::NumTUcells() {
	return TUcells.NumTUcells();
}

void SIMcore::getState(mxArray*** theOutput) {
	//this function outputs the inner data to MATLAB

	//1 - creating the most outer structure
	const char *field_names[] = { "env", "TUcells","Lymphocytes" }; 
	mwSize dims[2] = { 1, 1 };
	(*theOutput)[0] = mxCreateStructArray(1, dims, 3, field_names); //4

	//exporting environment
	mxSetField((*theOutput)[0],0,"env",env.getState());

	//exporting TUcells
	mxSetField((*theOutput)[0],0,"TUcells",TUcells.getState());

	//exporting TUcells
	mxSetField((*theOutput)[0], 0, "Lymphocytes", lymphocytes.getState());

}

void SIMcore::TUcellsNum(mxArray*** theOutput) {
	//this function outputs the inner data to MATLAB
	mwSize dims[3] = { 1, 1, 1 };
	//1 - creating the most outer structure
	(*theOutput)[0] =  mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
	int *ptr = (int*)mxGetPr((*theOutput)[0]);
	ptr[0] = (int)TUcells.NumTUcells();		
}
