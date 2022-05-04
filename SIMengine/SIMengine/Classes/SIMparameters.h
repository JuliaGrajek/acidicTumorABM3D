/*
 *  Created by Jan Poleszczuk
 *  Last modified 2022 by J Grajek
 */


#pragma once
#include "stdafx.h"

class SIMparameters
{
public:
    
    double ATPthresh; //cells die when they do not produce sufficient ATP, i.e. when nutrients are scarce
	//tumor cell related parameters
	unsigned char TUpmax; //divisions before proliferative capacity exhaustion
    float PDL1SuppProb;
    double keSLC0111;
    double kaSLC0111;
    double keantiPD1;
    double kaantiPD1;
    
    double CellDieAtProl; // cell dies at proliferation attempt, default 0, increase by chemo
	double TUpprol; //division probability
	double TUpdeath; //spontaneous death probability
	double TUps; //probability of symmetric division
	double TUpmig; //probability of migration
	double TUpmut; //probability of mutation
	float TUdanti; //antigenicity strength of mutating tumor cell
    double TUprotThresh; // proton level at which TC die with prob=1
    double TUARprotThresh; //proton level at which acidresistant TC die with prob=1
    double TUprotThreshQuiescence; // proton level at which TC undergo quiescence, i.e. survive but can't proliferate (reversible)

	float defaultAntigenicity;
	double stromaPerm;
	float maxAntigenicity;

	float IMIFNg; //adjuvanticity strength of dying tumor cell
	int IFNgRange;
	float IFNgDecay;

	int smoothRadius;
    bool seedUnderneath;
	
	//chemotaxis map associated parameters
	double DCchemo;
	double SCchemo;
    
    //proton map associated parameters
	double DCproton;
	double SCproton;
    double pHBuffer;
    double physiologicalProton;
    
	//necrosis map associated parameters
	double DCnecro;
	double TCnecro;
    double hypThresh;
    double physiologicalOxygen;
    double oxygenPointConsumption;
    
    //glucose map associated parameters
    double DCglucose;
    double TCglucose;
    double glucThresh;
    double GlycTumRate;
    double physiologicalGlucose;
    double glucosePointConsumption;

	//lymphocytes assiociated parameters
	unsigned char IMpmax;
	unsigned char IMkmax;
	double IMpmig;
	double IMpprol;
	double IMpdeath;
	float IMrwalk;
	int IMinfluxRate;
	int IMspeed;
	double probSeedFibr;
	double fibrFrac;
    double IMhypoDeath;
    double IMprotThresh;
    double IMprotThreshQuiescence;
    float CA9sup;
    float CA9freq;
    double CA9protons;
    float PDL1freq;

	float antiTresh;
	float IFNgThresh;

	unsigned char damageTresh;
	unsigned char engagementDuration;

	unsigned int N1, N2, N3;; //environment dimensions
    float IFNgMap;
    bool Lf;
    
	double initSeed;

	SIMparameters();
	~SIMparameters();

	void listSIMParameters();

	void initialize(const mxArray*, const mxArray*);

};

