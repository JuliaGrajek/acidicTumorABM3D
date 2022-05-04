% JN Kather, NCT Heidelberg, 2017
%updated by J Grajek 2022

function [mySystem, cnst] = getSystemParams(spaceSetting)
disp('requested system parameters');

% START GENERAL SYSTEM PROPERTIES -------------------------------------------
mySystem.params.initialSeed = [];        % initial random seed, default 1
mySystem.params.seedUnderneath = false; % if fibrosis can be seeded below cell
mySystem.params.CellDieAtProl = 0; % cell dies upon proliferation attempt, default 0 (increase by chemo)
mySystem.params.PDL1SuppProb=1; %probability that immune cell is surpressed when attacking a PDL1 expressing tumor cell, default 1 (decreases with anti PD1/PDl1 treatment)
% END SYSTEM PROPERTIES -------------------------------------------

% START INITIALIZE TUMOR CELLS -------------------------------------------
mySystem.params.TUpprol = 0.5055;   % HISTO GENERATED - probability of proliferation
mySystem.params.TUpmig = 0.35;      % probability of migrating, default 0.35 (<1-TUpprol-TUpdeath)
mySystem.params.TUpdeath =  1-(1-0.0319)^4;  % HISTO GENERATED - probability of spontaneous death
mySystem.params.TUpmax = 10;        % FIXED max. proliferation capacity, default 10
mySystem.params.TUdanti = 0.1;      % FIXED antigenicity strength of mutating tumor cell
mySystem.params.TUdamageThresh = 2; % T cell inflicted damage threshold for tumor cell, default 2
mySystem.params.TUps = 0.42;    % 3D probability of symmetric division, default 0.42
mySystem.params.TUpmut = 0.4;     % mutation probability (increases antigenicity)



% END INITIALIZE TUMOR CELLS ---------------------------------------------

% START INITIALIZE LYMPHOCYTES ------------------------------------------
mySystem.params.IMkmax = 10;          % FIXED killing capacity of immune cells, default 10
mySystem.params.IMpmax = 10;          % FIXED proliferation capacity of immune cells, default 10
mySystem.params.IMpmig = 0.3;       % probability of lymphocyte migration, default 0.3
mySystem.params.IMrwalk = 0.5;       % random influence on movement, default 0.5
mySystem.params.IMspeed = 97;        % speed of immune cell movement, default 97
mySystem.params.IMpprol = 0.0449/mySystem.params.IMspeed;   % HISTO GENERATED - probability of proliferation
mySystem.params.IMpdeath = (1-(1-0.0037)^4)/mySystem.params.IMspeed;  % HISTO GENERATED - probability of spontaneous death
mySystem.params.engagementDuration = 48; % how many intermed. steps is a killing cell engaged? default 48 (=6 hrs)
mySystem.params.IFNgThresh = 1e-4;     % FIXED IFNgamma threshold for PD-L1 induction, default 0.2 
mySystem.params.antiThresh =0.3;    % FIXED antigenicity threshold for lymphocyte activation, default 0.3 
mySystem.params.IFNgDecay = 0.1;     % IFNgamma decay in each iteration
mySystem.params.IFNgRange = 2;      	% effect range of TC on IFNgMAP (radius of circle), default 7
mySystem.params.IMIFNg = 1.02*1e-4; %IFNg produced by active T cells
mySystem.params.IMinfluxProb = 0.33; % probability of immune cell influx, default 0.33
mySystem.params.IMhypoDeath = 1.2; %increase in apoptosis probability under hypoxic conditions
mySystem.params.IMinfluxRate = 1; 	   % 3D how many lymphocytes appear simultaneously, always fixed at 1
mySystem.params.IMrateDynamic = 0.04;   % 3D how does lymphocyte influx scale with increasing tumor size
% END INITIALIZE LYMPHOCYTES --------------------------------------------


% START INITIALIZE CHEMOTAXIS MAP ------------------------------------------
mySystem.params.DCchemo = 100; %diffusion/consumption in the stationary diffusion-consumption equation
mySystem.params.SCchemo = 1; %secretion/consumption in the stationary diffusion-consumption equation
% END INITIALIZE CHEMOTAXIS MAP  ---------------------------------


% START INITIALIZE NECROSIS
oxygenDiffusion = 2.5*1e-5*12*60*60; %ref (Powathil, et al., Comput Math Met Med, 2012), in cm^2/12h
mySystem.params.oxygenPointConsumption = 3.8*1e-13*12*60*60; %ref (Powathil, et al., Comput Math Met Med, 2012), in cm^2*O2/(12h*cell)
carryingCapacity = 2.1*10^11/2; %ref (Powathil, et al., Comput Math Met Med, 2012); scaled by 2, cells
dx = 14.9*1e-4;%grid spot width, in cm
mySystem.params.DCnecro = oxygenDiffusion/mySystem.params.oxygenPointConsumption/carryingCapacity/(3*dx)^2; %diffusion/(TU cell consumption) in the stationary diffusion-consumption equation + numerical scheme correction
mySystem.params.TCnecro = 1; % FIXED total consumption by a single TU cell (because of normalization D/c)
mySystem.params.physiologicalOxygen = 0.056; %in mM, oxygen concentration in blood
mySystem.params.hypThresh = mySystem.params.physiologicalOxygen/12; %threshold of oxygen value, below which we assume hypoxia.
% END INITIALIZE NECROSIS

% START INITIALIZE GLUCOSE 
% parameter values
mySystem.params.physiologicalGlucose =  5; %Glucose in Blood
glucoseDiffusion = 2.6*1e-6*12*60*60; %glucose diffusion coefficient
mySystem.params.glucosePointConsumption = 10/3*mySystem.params.oxygenPointConsumption*mySystem.params.physiologicalOxygen/(5*mySystem.params.physiologicalGlucose); %assumption: in physiological conditions 70% of glucose is used for glycolysis TCs only uptake oxygen
%carryingCapacity = 2.1*10^11/2; %ref (Powathil, et al., Comput Math Met Med, 2012); scaled by 2, cells
mySystem.params.DCglucose =glucoseDiffusion/mySystem.params.glucosePointConsumption/carryingCapacity/(3*dx)^2; %diffucion/(TU cell consumption) in the stationary diffusion-consumption equation + numerical scheme correction
mySystem.params.TCglucose = 1; % FIXED total consumption by a single TU cell (because of normalization D/c)
mySystem.params.glucThresh = 0.5; % threshold of nutrients value below which lymphocytes have an impaired effector function (cytotoxicitiy) and motility
mySystem.params.GlycTumRate = 9/7; %glyc_rate for TC in hypoxic regions
% END INITIALIZE  GLUCOSE

% ATP
mySystem.params.ATPthresh =(mySystem.params.physiologicalOxygen*mySystem.params.oxygenPointConsumption*29/5)*0.05; %ATP threshold below which cells die. 

% PK
mySystem.params.keSLC0111 = log(2)/(11.1/12); %half-life of SLC-0111 is 11.1h, one time step is 12h
mySystem.params.kaSLC0111 = 0.93*12; %calculated using tmax and the above calculated ke
mySystem.params.keantiPD1 = log(2)/(26*2); %elimination half life is 26 days
mySystem.params.kaantiPD1 = 2.6*12; %not used remove
%END PK

% START INITIALIZE Proton MAP ------------------------------------------
mySystem.params.pHBuffer = 5*1e-5; %parameter encompassing all of the proton buffering
protonDiffusion =  1.08*1e-5*12*60*60; %proton diffusion coefficient
protonPointSecretion = 2*mySystem.params.pHBuffer; %coefficient of proton secretion. In ProtonMap.cpp it is additionally
%multiplied by the amount of glucose uptaken that wasn't used for aerobic respiration
mySystem.params.DCproton = protonDiffusion/protonPointSecretion/carryingCapacity/(3*dx)^2; %diffusion/consumption in the stationary diffusion-consumption equation, 
mySystem.params.SCproton = 1; %secretion/consumption in the stationary diffusion-consumption equation 
mySystem.params.physiologicalProton = 3.98*1e-5; %physiological H+ concentration in mM, corresponding to pH=7.4
mySystem.params.TUprotThresh = 1e-3; %proton threshold in mM at which TCs die. Equivalent to pH=6 
mySystem.params.TUARprotThresh = 1.5*1e-3; %proton threshold at which acidresistant TCs die Equivalent to pH=5.8.
mySystem.params.TUprotThreshQuiescence = 3.98*1e-4; %proton threshold at which TCs undergo quiescence. Equivalent to pH=6.4.
mySystem.params.IMprotThresh = 3.98*1e-4; %proton threshold in mM at which IMs die Equivalent to pH=6.4 
mySystem.params.IMprotThreshQuiescence = 2*1e-4; %proton threshold at which  IMs undergo quiescence. Equivalent to pH=6.7. 
mySystem.params.CA9sup = 0; %probability of inhibiting CAIX (positive when adding SLC-0111 treatment)
mySystem.params.CA9freq =0.3; %frequency of CAIX expression
mySystem.params.CA9protons =2.1*1e-9;%4*1e-9; %how many protons more are produced due to CAIX
mySystem.params.PDL1freq=0; %initial PDL1 frequency
% END INITIALIZE Proton MAP  ---------------------------------

% START INITIALIZE FIBROSIS  ---------------------------------
mySystem.params.smoothRadius = 3; 		% for smoothing fibrotic maps, default 3
mySystem.params.probSeedFibr = 0.06;   	% probability of turning into fibrosis, model fitting: 0.06
mySystem.params.fibrFrac = 0.3;    		% FIXED size of fibrotic seed, 0...1, default 0.3
mySystem.params.stromaPerm = 1;	% 0 = stroma not permeable, 1 = fully permeable
% END INITIALIZE FIBROSIS  ---------------------------------

% START DEFINING ADDITIONAL CONSTANTS -----------------------------------
cnst.verbose = true;            % draw intermediary steps? default true
cnst.createNewSystem = true;    % create new system at the start, default true
cnst.maxCells = Inf; %0.04*pi/(dx)^2; %maximal number of tumor cells
cnst.averageOut = true;

% cnst.saveImage = true;          % save image directly after each iteration, default true
% cnst.doImage = false;           % plot result again afterwards, default false
% cnst.doVideo = false;           % create a video afterwards, default false
% cnst.doSummary = true;          % summarize the result, default true
cnst.inTumor = 1;               % defines "in tumor" ROI, default 1
cnst.marginSize = round(67/2);  % default "invasive margin" ROI, default 67/2 = 0.5 mm
cnst.around = round(67*2);      % defines "adjacent tissue" ROI, default 67*2 = 2 mm
cnst.requireAlive = 150;        % require tumor to be alive for some time
cnst.maxAntigenicity = 1;       % maximum antigenicity for tumor cells
cnst.tumorColorLevels = 100;    % how many tumor color shades (must be multiple of 4!)
cnst.antigenKernelWidth = 0.05; % width of ksdensity kernel for antigenicity plot
cnst.defaultAntigenicity = 0.05; % antigenicity of first tumor cell
cnst.smoRegion = strel('disk',5,0); % region smoother for topography statistics
cnst.topoColorNorm = 0.03;      % max value for topography maps (red/blu)
cnst.VideoFrameRep = 5;         % repetition of video frames
cnst.lossFunction = 'default';  % define the loss function for system summary, 'default' or 'none'
cnst.penalty = 50;              % define penalty for immune cell win for loss function, default 50
% END DEFINING ADDITIONAL CONSTANTS -----------------------------------

% START DEFINING DIMENSION VARIABLES  -----------------------------------

mySystem.is3D = true;
mySystem.grid.N = spaceSetting(1);  % domain dimension 1, default 90
mySystem.grid.M = spaceSetting(2);  % domain dimension 2, default 90
mySystem.grid.P = spaceSetting(3);  % domain dimension 3, default 90
mySystem.grid.IFNgMap = zeros(mySystem.grid.N, mySystem.grid.M, mySystem.grid.P, 'single');
mySystem.grid.Lf = logical(zeros(mySystem.grid.N, mySystem.grid.M, mySystem.grid.P, 'single'));


% START PLAUSIBILITY CHECK ----------------------------------------------
checkPlausibility(mySystem);
% END PLAUSIBILITY CHECK ------------------------------------------------

end
