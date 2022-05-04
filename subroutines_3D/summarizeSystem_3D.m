% JN Kather 2017, jakob.kather@nct-heidelberg.de
% updated by J Grajek, 2022
%
% compatible with model 2.0(TU/IM/MP):     yes
% compatible with 3D:                      no

function summaryOut = summarizeSystem_3D(mySystem,cnst)
%summarizeSystem_NM summarizes the state of a 3D system
%   Input is the "mySystem" structure with fields grid, params, TU, IM
%   and the "cnst" structure, containing basic constants
%   This function genreates a structure "summaryOut" containing all 
%   relevant measurements of the system's state
    
    %Find probability that there is a T cell in the vicinity of a tumor
    %cell
    neigh=[];
    N1=mySystem.grid.N;
    N2=mySystem.grid.M;

    for i=-2:2
        for j=-2:2
            for k=-2:2
                neigh=[neigh, k*N1*N2+j*N1+i];
            end
        end
    end

    TILsVicinity=0;

    
    for c=1:numel(mySystem.TU.TUcells)
        TILsToCheck = double(mySystem.TU.TUcells(c))+neigh;
        if intersect(TILsToCheck, mySystem.IM.IMcells)
            TILsVicinity=TILsVicinity+1;
        end
    end
    
    % create basic results
    summaryOut.TU_Num = numel(mySystem.TU.TUcells); % tumor cell number
    summaryOut.TU_MaxDiameter = [];
    summaryOut.TU_FracStem = sum(mySystem.TU.TUprop.isStem) ...
        / numel(mySystem.TU.TUprop.isStem); % stem cell fraction
    summaryOut.IM_Num = sum(mySystem.IM.IMprop.Kcap > 0 & mySystem.IM.IMprop.quiescent==0); % immune cell number
    summaryOut.IM_FracExhaust = sum(mySystem.IM.IMprop.Kcap == 0); % fraction of exhausted immune cells
    summaryOut.IM_FracQuiescent = sum(mySystem.IM.IMprop.quiescent == 1);
    summaryOut.minpH = min(-log10(mySystem.grid.ProtMap/1000), [], 'all');
    summaryOut.minO2 = min(mySystem.grid.HypoxMap, [], 'all');
    summaryOut.minG = min(mySystem.grid.GlucMap, [], 'all');
    summaryOut.minATP = min(mySystem.grid.ATPMap, [], 'all');
    summaryOut.PDLfrac = sum(mySystem.TU.TUprop.PDL)/numel(mySystem.TU.TUprop.PDL);
    summaryOut.meanpH=mean(-log10(mySystem.grid.ProtMap(mySystem.TU.TUcells)/1000),'all');
    summaryOut.CA9freq = sum(mySystem.TU.TUprop.isAcidResistant)/numel(mySystem.TU.TUprop.isAcidResistant);
    summaryOut.PDL1SuppProb = mySystem.params.PDL1SuppProb;   

	% if there is a chemotaxis map, then create spatial results
    if isfield(mySystem.grid,'ChtaxMap')
    ChtaxMap2 = double(bwdist(mySystem.grid.Lt,'euclidean')); 
    Mask_fullTumorAux = imfill(ChtaxMap2<cnst.inTumor,'holes');   % binary mask ROI 1
    Mask_tumorCore = imerode(Mask_fullTumorAux,strel('disk',round(cnst.marginSize/2),6));
    Mask_marginIn = Mask_fullTumorAux & ~Mask_tumorCore;
    Mask_marginOut = ChtaxMap2<cnst.marginSize & ~Mask_fullTumorAux; % binary mask ROI 2
    Mask_distantOut = ChtaxMap2<cnst.around & ~Mask_marginOut & ~Mask_fullTumorAux;  % binary mask ROI 3
	
	% count lymphocyte (IM) cells in regions PER GRID CELL
    summaryOut.IM_tumorCore = sum(Mask_tumorCore(mySystem.IM.IMcells))/sum(Mask_tumorCore(:));
    summaryOut.IM_marginIn = sum(Mask_marginIn(mySystem.IM.IMcells))/sum(Mask_marginIn(:));
    summaryOut.IM_marginOut = sum(Mask_marginOut(mySystem.IM.IMcells))/sum(Mask_marginOut(:));
    summaryOut.IM_distantOut = sum(Mask_distantOut(mySystem.IM.IMcells))/sum(Mask_distantOut(:));
   
	% more features: tumor/stroma ratio, tumor/necrosis ratio, stromal lymphocytes fraction
    summaryOut.TU_Stro_ratio_log = log(double(numel(mySystem.TU.TUcells))/double(sum(mySystem.grid.Lf(:))));
    summaryOut.Stro_Fraction = double(sum(mySystem.grid.Lf(:)))/double(numel(mySystem.TU.TUcells));
    summaryOut.TU_Necr_ratio_log = log(double(numel(mySystem.TU.TUcells))/double(sum(mySystem.grid.Ln(:))));
    summaryOut.IM_instroma = sum(mySystem.grid.Lf(mySystem.IM.IMcells))/numel(mySystem.IM.IMcells);
    summaryOut.minpH = min(-log10(mySystem.grid.ProtMap/1000), [], 'all');
   
    
    % calculate immune escape indices: fraction of tumor cells in
    % low-antigenicity (sub-threshold) or low-adjuvanticity (sub-threshold)
    % conditions
    inflammatoryTME =  mySystem.grid.IFNgMap(mySystem.TU.TUcells); % get all IFNgamma values in Tumor Cell locations
    summaryOut.immuneEscape.lowAntigenicity  = sum(mySystem.TU.TUprop.Antigen <= mySystem.params.antiThresh) / numel(mySystem.TU.TUcells);
    summaryOut.immuneEscape.meanAntigenicity = mean(mySystem.TU.TUprop.Antigen);
    summaryOut.immuneEscape.lowIFNg = sum( inflammatoryTME <= mySystem.params.IFNgThresh)  / numel(mySystem.TU.TUcells); 
    
        if cnst.verbose
            disp(['immune escape measurements: ',10,'low anti: ',num2str(round(100*summaryOut.immuneEscape.lowAntigenicity)),'%',...
                ', low adju: ', num2str(round(100*summaryOut.immuneEscape.lowIFNg)),'%']);
        end
    end
    
    % copy hyper-parameters
    summaryOut.stepsDone = mySystem.grid.StepsDone;
    summaryOut.N = mySystem.grid.N;
    summaryOut.M = mySystem.grid.M;
    
    summaryOut.TILsVicinityProb = TILsVicinity/summaryOut.TU_Num;
    
end