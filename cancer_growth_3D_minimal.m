
% 2016-2022, created by JN Kather and J Poleszczuk and J Grajek
% license: see separate license file

% What is this? this is an agent-based model of tumor cell - immune cells
% interactions, which takes into account interactions with the
% tumor-microenvironment (TME).

% How does it work? open this file in Matlab and run it.

% the code relies on a simulation engine that it written in C++. It has
% been compiled on a Windows computer. Before you run the model on a 
% Linux or MacOS computer, you might have to re-compile it. 
% See instructions in ./SIMengine
% 
% the C++ code relies on the open source software SuiteSparse (DLLs are
% included in the ./SIMengine folder. For license and disclaimer please see
% https://github.com/jlblancoc/suitesparse-metis-for-windows/

close all; clear variables; format compact; clc % avoid spillover 
addpath('./SIMengine/'); % include SIMengine (MEX-based simulation engine)
addpath('./subroutines/'); % include generic subroutines for 2D and 3D
addpath('./subroutines_3D/'); % include generic subroutines for 3D modeling
addpath('./subroutines_plot/'); % include advanced subroutines for plotting

numExp =1;                 % number of simulation runs, default 1
saveImage = true;           % save simulation output image. requires verbose true
saveVideo = true;          % save simulation output video. requires verbose true

%Radius of initial tumor, if you don't want to start simulation from a
%single cell, but from a sphere in the center of the domain
r=30; %default 30 cells

% all parameters for the model are stored in the structure "sysTempl".
% Hyperparameters are stored in the structure "cnst". If you want to 
% manually change parameters, you need to overwrite the respective value 
% in sysTempl.params or in cnst, for example by adding
% "sysTempl.params.TUps = 0.65" after the call to "getSystemParams"

for i=1:numExp
    [sysTempl, cnst] = getSystemParams([180 180 180]);  % get system parameters. 2nd argument is domain size which must be
                                                             % multiplication of 3. [135 135 135] is a 2 mm cube domain.

    cnst.nSteps   = 80; %how many iterations. 1 iteration = 12 hours
    cnst.drawWhen = 5;  % update plot after ... iterations
    
    %for pharmacokinetics (PK) simulations
    cnst.nStepsBeforeTreatment = 10; %when simulating PK, how long does the tumor grow before administering anti-PD1
    cnst.nStepsAfterTreatment = 14*2; %when simulating PK, how long does the tumor grow after administering anti-PD1

    % add/override some global variables after loading the system

    cnst.VideoFrameRep = 12;    % frame repetition if video is recorded, default 12
    cnst.verbose = true;        % show simulation output on screen
    cnst.printImages = true;   % save high resolution simulation output image
    cnst.lossFunction = 'default_stem'; % specify loss function for simulation
    sysTempl.experiment_id=2;            %1: in vivo from scratch, 2: in vivo from big sphere of radius r, 3: with pharmakodynamics
    sysTempl.lym_id=0;
    sysTempl.params.initialSeed = i; % reset random seed for reproducibility
    expname = ['minimal_3D_',num2str(i)]; % experiment name for saving
    
    %SIMULATION
    %sysTempl.params.PDL1SuppProb=0.75;
    [sysOut, lastFrame, summary, imWin, masterID] = ...
        runSystem(sysTempl,cnst,expname,saveImage,saveVideo, sysTempl.experiment_id, r);
    
    
end

