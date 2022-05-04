/*
*  Created by Jan Poleszczuk
*  Last modified 2022 by J Grajek
*/


#include "SIMparameters.h"


SIMparameters::SIMparameters()
{
	// START GENERAL SYSTEM PROPERTIES -------------------------------------------
	initSeed = 0;
	seedUnderneath = true;
    CellDieAtProl = 0;
    PDL1SuppProb = 1;
	// END SYSTEM PROPERTIES -------------------------------------------
    ATPthresh = 0.5;
    keSLC0111 = log(2)/(9.8/12);
    kaSLC0111 =1.24*12;
    keantiPD1 = log(2)/(26.7*2);
    kaantiPD1 =2.6*12;
	// START INITIALIZE TUMOR CELLS -------------------------------------------
	TUpmax = 10; //divisions before proliferative capacity exhaustion
	TUpprol = 1. / 24.; //division probability
	TUpdeath = 0.05; //spontaneous death probability
	TUps = 0.3; //probability of symmetric division
	TUpmig = 10. / 24.; //probability of migration
	TUpmut = 0.;
	TUdanti = 0.1f;
    damageTresh = 5;
    TUprotThresh = 1;
    TUARprotThresh = 1;
    TUprotThreshQuiescence = 2;
    
	// END INITIALIZE TUMOR CELLS ---------------------------------------------


	// START INITIALIZE LYMPHOCYTES ------------------------------------------
	IMkmax = 5;
	IMpmax = 10;
	IMpmig = 0.8;
    IMIFNg = 1.f;
	IMrwalk = 0.8f;
	IMspeed = 30;
	IMpprol = 0.049 / (double)IMspeed;
	IMpdeath = 0.0147;
	engagementDuration = 48;
	IFNgThresh = 0.0f;
	antiTresh = 0.3f;
	IFNgDecay = 0.95f;
	IMinfluxRate = 1;
    IMhypoDeath = 1.2; 
    IMprotThresh = 1;
    IMprotThreshQuiescence =2;
	// END INITIALIZE LYMPHOCYTES --------------------------------------------

	// START INITIALIZE CHEMOTAXIS MAP------------------------------------------
	DCchemo = 1440.;
	SCchemo = 10.;
	// END INITIALIZE CHEMOTAXIS MAP------------------------------------------
       
	// START INITIALIZE PROTON MAP------------------------------------------
	DCproton = 1440.;
	SCproton = 10.;
    pHBuffer = 0.8;
    physiologicalProton = 1.;
    CA9freq = 0.09;
    CA9sup=0;
    CA9protons = 1;
    PDL1freq=0;
	// END INITIALIZE PROTON MAP------------------------------------------
	
	// START INITIALIZE NECROSIS------------------------------------------
	DCnecro = 1.;
	TCnecro = 1.;
    hypThresh = 0.5;
    physiologicalOxygen = 1.;
    oxygenPointConsumption =1.;
	// END INITIALIZE NECROSIS------------------------------------------
    
    //START INITIALIZE GLUCOSE
    DCglucose = 1.;
    TCglucose = 1.;
    glucThresh = 0.5;
    GlycTumRate = 5.;
    physiologicalGlucose = 1.;
    glucosePointConsumption = 1.;
	//END INITIALIZE GLUCOSE
    
	// START INITIALIZE FIBROSIS  ---------------------------------
	smoothRadius = 3;
	probSeedFibr = 0.008;
	fibrFrac = 0.3;
	stromaPerm = 0.1;
	// END INITIALIZE FIBROSIS  ---------------------------------
	
	defaultAntigenicity = 0.05f;
	maxAntigenicity = 1.;
	

	N1 = 351; N2 = 351; N3 = 1;
//     A = 0.1f; F= 0;
}


SIMparameters::~SIMparameters(){};


void SIMparameters::listSIMParameters() {
	printf("Current simulation parameters in SIMengine: \n");
	printf("_____________________________________________ \n");
	printf("TUpmax: %d \n", TUpmax);
	printf("TUpprol: %f \n", TUpprol);
	printf("TUpdeath: %f \n", TUpdeath);
	printf("TUps: %f \n", TUps);
	printf("TUpmig: %f \n", TUpmig);
	printf("TUpmut: %f \n", TUpmut);
	printf("TUdanti: %f \n", TUdanti);
	printf("IMIFNg: %f \n", IMIFNg);
	printf("IFNgRange: %d \n", IFNgRange);
	printf("IFNgDecay: %f \n", IFNgDecay);
	printf("DCnecro: %f \n", DCnecro);
	printf("DCchemo: %f \n", DCchemo);
	printf("SCchemo: %f \n", SCchemo);
	printf("defaultAntigenicity: %f \n", defaultAntigenicity);
	printf("maxAntigenicity: %f \n", maxAntigenicity);
	printf("IMpmax: %d \n", IMpmax);
	printf("IMkmax: %d \n", IMkmax);
	printf("IMpmig: %f \n", IMpmig);
	printf("IMpprol: %f \n", IMpprol);
	printf("IMpdeath: %f \n", IMpdeath);
	printf("IMrwalk: %f \n", IMrwalk);
	printf("IMinfluxRate: %d \n", IMinfluxRate);
	printf("IMspeed: %d \n", IMspeed);
	printf("probSeedFibr: %f \n", probSeedFibr);
	printf("fibrFrac: %f \n", fibrFrac);
	printf("damageTresh: %d \n", damageTresh);
	printf("engagementDuration: %d \n", engagementDuration);
	printf("antiTresh: %f \n", antiTresh);
	printf("IFNgThresh: %f \n", IFNgThresh);
// 	printf("MPpmax: %d \n", MPpmax);
// 	printf("MPinfluxRate: %d \n", MPinfluxRate);
	printf("stromaPerm: %f \n", stromaPerm);
	//printf("initSeed: %f \n", initSeed);
	printf("smoothRadius: %d \n", smoothRadius);
	printf("N1xN2xN3: %dx%dx%d \n", N1, N2, N3);
}

void SIMparameters::initialize(const mxArray* sysTempl, const mxArray* cnst) {
	//parsing input parameters
	mxArray *params = mxGetField(sysTempl, 0, "params");
    ATPthresh= (double)mxGetPr(mxGetField(params, 0, "ATPthresh"))[0];
    CellDieAtProl = mxGetPr(mxGetField(params, 0, "CellDieAtProl"))[0];
    PDL1SuppProb = mxGetPr(mxGetField(params, 0, "PDL1SuppProb"))[0];
	TUpmax = (unsigned char)mxGetPr(mxGetField(params, 0, "TUpmax"))[0]; //divisions before proliferative capacity exhaustion
	TUpprol = mxGetPr(mxGetField(params, 0, "TUpprol"))[0]; //division probability
	TUpdeath = mxGetPr(mxGetField(params, 0, "TUpdeath"))[0]; //spontaneous death probability
	TUps = mxGetPr(mxGetField(params, 0, "TUps"))[0]; //probability of symmetric division
	TUpmig = mxGetPr(mxGetField(params, 0, "TUpmig"))[0]; //probability of migration
	TUpmut = mxGetPr(mxGetField(params, 0, "TUpmut"))[0]; //probability of mutation
	TUdanti = (float)mxGetPr(mxGetField(params, 0, "TUdanti"))[0];
	IMIFNg = (float)mxGetPr(mxGetField(params, 0, "IMIFNg"))[0];
    TUprotThresh = (double)mxGetPr(mxGetField(params, 0, "TUprotThresh"))[0];
    TUARprotThresh = (double)mxGetPr(mxGetField(params, 0, "TUARprotThresh"))[0];
    TUprotThreshQuiescence = (double)mxGetPr(mxGetField(params, 0, "TUprotThreshQuiescence"))[0];
    keSLC0111 = (double)mxGetPr(mxGetField(params, 0, "keSLC0111"))[0];
    kaSLC0111 = (double)mxGetPr(mxGetField(params, 0, "kaSLC0111"))[0];
    keantiPD1 = (double)mxGetPr(mxGetField(params, 0, "keantiPD1"))[0];
    kaantiPD1 = (double)mxGetPr(mxGetField(params, 0, "kaantiPD1"))[0];
    
	initSeed = mxGetPr(mxGetField(params, 0, "initialSeed"))[0];
	stromaPerm = (double)mxGetPr(mxGetField(params, 0, "stromaPerm"))[0];

	IFNgRange = (int)mxGetPr(mxGetField(params, 0, "IFNgRange"))[0];
	IFNgDecay = (float)mxGetPr(mxGetField(params, 0, "IFNgDecay"))[0];

	smoothRadius = (int)mxGetPr(mxGetField(params, 0, "smoothRadius"))[0];
    seedUnderneath = mxGetPr(mxGetField(params, 0, "seedUnderneath"))[0] > 0;

	DCnecro = mxGetPr(mxGetField(params, 0, "DCnecro"))[0];
	TCnecro = mxGetPr(mxGetField(params, 0, "TCnecro"))[0];
    hypThresh = mxGetPr(mxGetField(params, 0, "hypThresh"))[0];
    physiologicalOxygen = mxGetPr(mxGetField(params, 0, "physiologicalOxygen"))[0];
    oxygenPointConsumption = mxGetPr(mxGetField(params, 0, "oxygenPointConsumption"))[0];
    
	DCproton = mxGetPr(mxGetField(params, 0, "DCproton"))[0];
	SCproton = mxGetPr(mxGetField(params, 0, "SCproton"))[0];
    pHBuffer = mxGetPr(mxGetField(params, 0, "pHBuffer"))[0];
    physiologicalProton = mxGetPr(mxGetField(params, 0, "physiologicalProton"))[0];
    CA9freq = mxGetPr(mxGetField(params, 0, "CA9freq"))[0];
    CA9sup = mxGetPr(mxGetField(params, 0, "CA9sup"))[0];
    CA9protons = mxGetPr(mxGetField(params, 0, "CA9protons"))[0];
    PDL1freq = mxGetPr(mxGetField(params, 0, "PDL1freq"))[0];
    
	DCchemo = mxGetPr(mxGetField(params, 0, "DCchemo"))[0];
	SCchemo = mxGetPr(mxGetField(params, 0, "SCchemo"))[0];
    
    DCglucose = mxGetPr(mxGetField(params, 0, "DCglucose"))[0];
	TCglucose = mxGetPr(mxGetField(params, 0, "TCglucose"))[0];
	glucThresh = mxGetPr(mxGetField(params, 0, "glucThresh"))[0];
    GlycTumRate = mxGetPr(mxGetField(params, 0, "GlycTumRate"))[0];
    physiologicalGlucose = mxGetPr(mxGetField(params, 0, "physiologicalGlucose"))[0];
    glucosePointConsumption = mxGetPr(mxGetField(params, 0, "glucosePointConsumption"))[0];
    
	IMpmax = (unsigned char)mxGetPr(mxGetField(params, 0, "IMpmax"))[0];
	IMkmax = (unsigned char)mxGetPr(mxGetField(params, 0, "IMkmax"))[0];
	IMpmig = mxGetPr(mxGetField(params, 0, "IMpmig"))[0];
	IMpprol = mxGetPr(mxGetField(params, 0, "IMpprol"))[0];
	IMpdeath = mxGetPr(mxGetField(params, 0, "IMpdeath"))[0];
	IMrwalk = (float)mxGetPr(mxGetField(params, 0, "IMrwalk"))[0];
	probSeedFibr = mxGetPr(mxGetField(params, 0, "probSeedFibr"))[0];
	fibrFrac = mxGetPr(mxGetField(params, 0, "fibrFrac"))[0];
    IMhypoDeath = mxGetPr(mxGetField(params, 0, "IMhypoDeath"))[0];;
    IMprotThresh = (double)mxGetPr(mxGetField(params, 0, "IMprotThresh"))[0];
    IMprotThreshQuiescence = (double)mxGetPr(mxGetField(params, 0, "IMprotThreshQuiescence"))[0];
    
	IMinfluxRate = (int)mxGetPr(mxGetField(params, 0, "IMinfluxRate"))[0];
	IMspeed = (int)mxGetPr(mxGetField(params, 0, "IMspeed"))[0];

	damageTresh = (unsigned char)mxGetPr(mxGetField(params, 0, "TUdamageThresh"))[0];
	engagementDuration = (unsigned char)mxGetPr(mxGetField(params, 0, "engagementDuration"))[0];

	antiTresh = (float)mxGetPr(mxGetField(params, 0, "antiThresh"))[0];
	IFNgThresh = (float)mxGetPr(mxGetField(params, 0, "IFNgThresh"))[0];

	mxArray *grid = mxGetField(sysTempl, 0, "grid");
	N1 = (unsigned int)mxGetPr(mxGetField(grid, 0, "N"))[0];
	N2 = (unsigned int)mxGetPr(mxGetField(grid, 0, "M"))[0];
	if (mxGetField(grid, 0, "P") != NULL) {
		N3 = (unsigned int)mxGetPr(mxGetField(grid, 0, "P"))[0];
	}
	else {
		N3 = 1;
	}
    IFNgMap = (unsigned int)mxGetPr(mxGetField(grid, 0, "IFNgMap"))[0];
    Lf = (unsigned int)mxGetPr(mxGetField(grid, 0, "Lf"))[0];
}
