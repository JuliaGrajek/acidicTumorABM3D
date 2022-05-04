/*
*  Created by Jan Poleszczuk
*  Last modified May, 2022 by JGrajek
*/

#include "stdafx.h"
#include "TumorCells.h"
#include <math.h>

size_t TumorCells::NumTUcells() {
	return cells.size();
}

void TumorCells::initialize(SIMparameters *parSet, Environment *envSet, CRandomMersenne* rG) {
	env = envSet;
	randGen = rG;
	params = parSet;

	//intializing with single cancer stem cell in the lattice center
	unsigned int center = env->getCenter();
	env->newTUcell(center, 1);
	TumorCell initialCell = { center, params->TUpmax, true, params->defaultAntigenicity, 0, 1, 0, 0 };
	cells.clear(); //just in case if simulation is reinitialized
	cells.push_back(initialCell);
	
}

void TumorCells::initializeFromState(const mxArray* sysTempl, SIMparameters *parSet, Environment *envSet, CRandomMersenne* rG) {
	env = envSet;
	randGen = rG;
	params = parSet;

	cells.clear(); //just in case

	//reading from state
	mxArray *TU = mxGetField(sysTempl, 0, "TU");

	mxArray *val = mxGetField(TU, 0, "TUcells");

	const mwSize *dimensionPtr = mxGetDimensions(val);
	int numTUCells = (int)dimensionPtr[1]; //assuming it is column vector

	unsigned int* pos = (unsigned int*)mxGetPr(val);
	
	mxArray *TUprops = mxGetField(TU, 0, "TUprop");
	val = mxGetField(TUprops, 0, "isStem");
	bool* iS = (bool*)mxGetPr(val);
	val = mxGetField(TUprops, 0, "Pcap");
	unsigned char* Pc = (unsigned char*)mxGetPr(val);
	val = mxGetField(TUprops, 0, "Antigen");
	float* An = (float*)mxGetPr(val);
	val = mxGetField(TUprops, 0, "damage");
	unsigned char* dem = (unsigned char*)mxGetPr(val);
    val = mxGetField(TUprops, 0, "glyc_rate");
	double* gr = (double*)mxGetPr(val);
    val = mxGetField(TUprops, 0, "isAcidResistant");
	bool* iA = (bool*)mxGetPr(val);
    val = mxGetField(TUprops, 0, "PDL");
	bool* PDL = (bool*)mxGetPr(val);

	for (int i = 0; i < numTUCells; ++i) {
		TumorCell cellTuAdd = { pos[i]-1, Pc[i], iS[i], An[i], dem[i], gr[i], iA[i], PDL[i] };
		env->newTUcell(pos[i]-1, gr[i]);
		cells.push_back(cellTuAdd);
	}
}

void TumorCells::readState(mxArray* init) {
	

}

void TumorCells::getStateForImmune(std::map<unsigned int, std::vector<TumorCells::TumorCell>::iterator>* TUmap) {

	std::vector<TumorCell>::iterator it;
	for (it = cells.begin(); it != cells.end(); it++) {
		(*TUmap)[it->place] = it;
	}

};

double TumorCells::metabolize(TumorCell currCell) { //TC metabolism: ATP and proton production based on the nutrient availability in the cell's location
            double gluc_uptake;
            double oxygen_uptake;
            double aerobic_glycolysis_glucose; //how much glucose remains after OXPHOS
            double ATP;
            
            gluc_uptake = env->getGlucoseValue(currCell.place)*currCell.glyc_rate*params->glucosePointConsumption; 
            oxygen_uptake = env->getNecrosisValue(currCell.place)*params->oxygenPointConsumption;
            
            //We assume that cells use OXPHOS as long as they have the oxygen to do it. OXPHOS uses 5 times as much oxygen as glucose. The remaining consumed glucose
            //is used for aerobic glycolysis. Glycolytic cells consume more glucose, so they use more aerobic glycolysis.
            aerobic_glycolysis_glucose = std::max(gluc_uptake-oxygen_uptake/5, static_cast<double>(0)); //net equation for aerobic glycolysis
            ATP = 2*aerobic_glycolysis_glucose + 29*std::min(oxygen_uptake/5, gluc_uptake);
            //from here aerobic glycolysis glucose describes the protons produced, so we add CA9protons 
            if (currCell.isAcidResistant){
                aerobic_glycolysis_glucose = aerobic_glycolysis_glucose+params->CA9protons/2; //divided by two, because glycolysis yields 2 protons and I included that in SCProton
            }
            //proton production: we assume the only source is aerobic glycolysis and CA9
            if (aerobic_glycolysis_glucose >0 && ATP>params->ATPthresh){
            env->newProtonSource(currCell.place, aerobic_glycolysis_glucose);
            }
            
            return ATP;

}


void TumorCells::go_grow_die() { 
    int killedCellsRound=0;
    unsigned char newCells=0;
	unsigned int newSite;
	TumorCell currCell, newCell;
	std::vector<TumorCell> cellsTmp;
	cellsTmp.reserve(cells.size());

	std::random_shuffle(cells.begin(), cells.end(), std::bind(&CRandomMersenne::Rand, *randGen, std::placeholders::_1)); //shuffling cells
	while (!cells.empty()) {
		currCell = cells.back(); //pick the cell
		cells.pop_back();
        
        if ((currCell.isAcidResistant==1) && (params->CA9sup!=0) &&(randGen->Random() < params->CA9sup)){

            currCell.isAcidResistant =0;            
        }
        
        
		if (currCell.damage == params->damageTresh || (!currCell.is_stem && randGen->Random() < params->TUpdeath)) {//cell is killed by the immune system, dies spontaneously
            env->deleteTUcell(currCell.place, currCell.glyc_rate);
            if (currCell.damage == params->damageTresh){
            killedCellsRound++;
            }
		}
        else if (env->getProtonValue(currCell.place)>params->TUprotThreshQuiescence && randGen->Random()<((env->getProtonValue(currCell.place)-params->TUprotThreshQuiescence)/((1-currCell.isAcidResistant)*params->TUprotThresh+currCell.isAcidResistant*params->TUARprotThresh-params->TUprotThreshQuiescence ))){ //cell dies due to low pH
            env->deleteTUcell(currCell.place, currCell.glyc_rate);
            env->markNecrosis(currCell.place);
        }           

        else{
            double ATP;
            ATP = metabolize(currCell);
            if (ATP < params->ATPthresh ){
                env->deleteTUcell(currCell.place, currCell.glyc_rate);
                env->markNecrosis(currCell.place);
        } 
            else{
            
            //Check for Hypoxia and adjust glycolytic rate accordingly
            if (env->getNecrosisValue(currCell.place) < params->hypThresh) {
                newCell = currCell; //we have to replace the cell because the glucose sink might have a different value if the cell was in a normoxic region when appearing and is now in a hypoxic region
                env->deleteTUcell(currCell.place, currCell.glyc_rate);
                newCell.glyc_rate = params->GlycTumRate;
                env->newTUcell(newCell.place, newCell.glyc_rate);
            }
            else if (env->getNecrosisValue(currCell.place) >= params->hypThresh){
                newCell = currCell; //we have to replace the cell because the glucose sink might have a different value if the cell was in a normoxic region when appearing and is now in a hypoxic region
                env->deleteTUcell(currCell.place, currCell.glyc_rate);
                newCell.glyc_rate = 1;
                env->newTUcell(newCell.place, newCell.glyc_rate);
           }
            if (env->getIFNgValue(currCell.place) > params->IFNgThresh){
                currCell.PDL=1;
            }
                    
			newSite = env->returnEmptyPlace(currCell.place);

			if (newSite) {//if there is a new spot
                newCell = currCell;
                newCell.place = newSite;
                newCell.PDL=0;
                
                if (randGen->Random() < params->TUpprol && env->getProtonValue(currCell.place)<params->TUprotThreshQuiescence) {
                    // START Chemotherapy
                    // tumor cell wants to divide
                    if (params->CellDieAtProl && (randGen->Random() < params->CellDieAtProl)) {
                        // tumor cell will die
                        env->deleteTUcell(currCell.place, currCell.glyc_rate);
                    }
                    // END Chemotherapy
                    else
                        {
                        if (currCell.is_stem) {
                            if (randGen->Random()<params->CA9freq) { //for now the CA9 distribution is random
                                newCell.isAcidResistant=1;
                            }
                            else{
                                newCell.isAcidResistant=0;
                            }
                            env->newTUcell(newSite, currCell.glyc_rate); //later maybe with a certain probability this, and o.w. 1
                            if (randGen->Random() > params->TUps) {//asymmetric division
                                newCell.is_stem = false;
                            }
                            if (currCell.Antigen < params->maxAntigenicity && randGen->Random() < params->TUpmut) {
                                currCell.Antigen += params->TUdanti;
                                if (currCell.Antigen > params->maxAntigenicity)
                                    currCell.Antigen = params->maxAntigenicity;
                                newCell.Antigen = currCell.Antigen;
                            }
                            
                            newCell.damage = 0;
                            cellsTmp.push_back(currCell);
                            cellsTmp.push_back(newCell);
                            newCells++;
                        }
                        else if (currCell.p==0){
                            env->deleteTUcell(currCell.place, currCell.glyc_rate);
                            
                        }
                        else{
                            //proliferation capacity is positive
                            if (randGen->Random()<params->CA9freq) { //for now the CA9 distribution is random
                                newCell.isAcidResistant=1;
                            }
                            else{
                                newCell.isAcidResistant=0;
                            }
                            currCell.p--;
                            newCell.p--;
                            env->newTUcell(newSite, currCell.glyc_rate); //later maybe inherits this with a cetain probabilit, ow. 1
                            if (currCell.Antigen < params->maxAntigenicity && randGen->Random() < params->TUpmut) {
                                currCell.Antigen += params->TUdanti;
                                if (currCell.Antigen > params->maxAntigenicity)
                                    currCell.Antigen = params->maxAntigenicity;
                                newCell.Antigen = currCell.Antigen;
                            }
                            newCell.damage = 0;
                            cellsTmp.push_back(currCell);
                            cellsTmp.push_back(newCell);
                            newCells++;
                            
                        }
                    }
                }
                else if (randGen->Random() < params->TUpmig) {
                    env->deleteTUcell(currCell.place, currCell.glyc_rate);
                    env->newTUcell(newSite, currCell.glyc_rate);
                    cellsTmp.push_back(newCell);
				}
				else {//doing nothing
                    
					cellsTmp.push_back(currCell);
				}
			}
			else {//no free spot
				cellsTmp.push_back(currCell);
			}
            }
		}
	}
	cells.swap(cellsTmp);
    
    killedCells.push_back(killedCellsRound);
    
}

mxArray* TumorCells::getState() {

	const char *field_names[] = { "position", "isStem", "Pcap", "Antigen","damage", "glyc_rate", "isAcidResistant", "PDL", "killed"
 };
	mwSize dimsStruct[2] = { 1, 1 };
	mxArray* TUcells = mxCreateStructArray(1, dimsStruct, 9, field_names);

	//2 - copying TUcells positions
    mwSize dims2[2];
    dims2[0]=1; dims2[1]=killedCells.size();
	mwSize dims[2];
	dims[0] = 1; dims[1] = cells.size();
	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "position"),
		mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL));
	int *TUcellsPos = (int*)mxGetPr(mxGetField(TUcells, 0, "position"));
	for (size_t i = 0; i < cells.size(); ++i)
		TUcellsPos[i] = cells.at(i).place+1;
    

	//3 - copying TUprop
	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "isStem"),
		mxCreateLogicalArray(2, dims));
	bool *isStem = (bool*)mxGetPr(mxGetField(TUcells, 0, "isStem"));

	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "Pcap"),
		mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL));
	unsigned char *Pcap = (unsigned char*)mxGetPr(mxGetField(TUcells, 0, "Pcap"));

	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "Antigen"),
		mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL));
	float *Antigen = (float*)mxGetPr(mxGetField(TUcells, 0, "Antigen"));

	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "damage"),
		mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL));
	unsigned char *dam = (unsigned char*)mxGetPr(mxGetField(TUcells, 0, "damage"));
    
    mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "glyc_rate"),
		mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL));
	double *glyc_rate = (double*)mxGetPr(mxGetField(TUcells, 0, "glyc_rate"));
    
    mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "isAcidResistant"),
		mxCreateLogicalArray(2, dims));
	bool *isAcidResistant = (bool*)mxGetPr(mxGetField(TUcells, 0, "isAcidResistant"));
    
    mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "PDL"),
    mxCreateLogicalArray(2, dims));
	bool *PDL = (bool*)mxGetPr(mxGetField(TUcells, 0, "PDL"));
    
    mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "killed"), mxCreateNumericArray(2, dims2,mxINT32_CLASS, mxREAL ));
    int *killed =(int*)mxGetPr(mxGetField(TUcells, 0, "killed"));
    for (size_t i=0; i < killedCells.size(); ++i){
        killed[i]=killedCells.at(i);
    }
    
	for (size_t i = 0; i < cells.size(); ++i) {
		isStem[i] = cells.at(i).is_stem;
		Pcap[i] = cells.at(i).p;
		Antigen[i] = cells.at(i).Antigen;
		dam[i] = cells.at(i).damage;
        glyc_rate[i] = cells.at(i).glyc_rate;
        isAcidResistant[i] = cells.at(i).isAcidResistant;
        PDL[i] = cells.at(i).PDL;
	}
    
	return TUcells;
}
