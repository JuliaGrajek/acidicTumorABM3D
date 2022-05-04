function [sysOut, lastFrame, summary, imWin, masterID] = ...
    runSystem(sysTempl,cnst,expname,saveImage,saveVideo, experiment_id, varargin)
% last preparations
masterID = ['Experiment_',char(expname)];
if (saveImage || saveVideo) && ~(exist(['./output/',masterID],'dir'))
mkdir(['./output/',masterID]); % create output directory to save results
end

%% RUN THE MODEL
% the model run is started by calling "growTumor_3D(...)". Because the
% tumor should not spontaneously die in the first few rounds, this process
% can be repeated four times until the tumor is still alive in the round
% specified by "requireAlive"
globalTime = tic;    % start timer
for nAttempts = 1:5  % try again if the tumor died too early

    if ~isnan(sysTempl.params.initialSeed)
        rng(sysTempl.params.initialSeed); % set random seed
    else
        rng('shuffle');
    end
    if experiment_id==1
        [sysOut, lastFrame, summary, imWin] = growTumor(sysTempl,cnst, experiment_id);
    elseif experiment_id==2
        [sysOut, lastFrame, summary, imWin] = growTumor(sysTempl,cnst, experiment_id, varargin{1});
    elseif experiment_id==3
        [sysOut, lastFrame, summary, imWin] = growTumorWithPK(sysTempl,cnst, experiment_id, varargin{1});
    else
        disp("Non-existent experiment id")
    end
            
% input variables:
% sysTempl: contains all parameters for the system (in the field
% sysTempl.params). Also, this strucutre will be used to save all model
% details later on. 
% cnst: contains all hyperparamters
%
% output variables:
% sysOut: this is a modified version of sysTempl that contains all
% information about the model in its final state, e.g. the positions and
% properties of all agents.
% lastFrame: contains a snapshot image of the the system everytime it was
% visualized (controlled by cnst.drawWhen)
% summary: contains summary statistics, e.g. number of tumor cells
% imWin: true if no tumor cells are left anymore
% fcount: how many times was the system visualized?

if imWin~=0 && imWin < cnst.requireAlive % attempt unsuccessful, try again with another seed
    initialSeed =initialSeed + 10;
else, break     % attempt successful
end

end
if cnst.verbose
    disp(['total time was ',num2str(toc(globalTime))]);
end

if cnst.printImages % print now in high quality - this causes parfor to crash
    print(gcf,['./output/p_',expname,'_',num2str(round(rand()*1000)),'.png'],'-dpng','-r450'); 
end

if saveImage % save last frame as image, low quality
    imwrite(lastFrame{end},['./output/',masterID,'/lastState_',num2str(sysTempl.params.initialSeed),'.png']); 
end
if saveVideo % save whole process as video, low quality for parfor, high quality for for
    writeMyVideo(lastFrame,['./output/',masterID,'/Movie_',num2str(sysTempl.params.initialSeed)],cnst.VideoFrameRep); 
end

end