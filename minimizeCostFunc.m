function [ out ] = minimizeCostFunc( in )
%minimizeCostFunc Uses simulated annealing to minimize STARE cost function
% See Bartlett et al. 2022 in NeuroImage for full STARE derivation to get
% equations that are contained in the constFunc.m and 2stModel.m that are
% called in this function.
%   Requires a structure as input with the following fields:
%   regions                             : Brain region names as strings in a cell array that will be quantified in STARE.
%   subOutPath                          : subject-level output directory. Created in stareAnchoring.m
%   vascCorrPerc                        : Percentage from 0 to 100 (0 and 5 used in validaions) for contribution to voxel-wise PET signal from vasculature signal (i.e., if =5, 5% of signal is corrupted from vasclature. For vascular correction of TACs.
%   midtime                             : string of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   weights                             : Generated frame-wise weights (square root of frame duration). Must be same size as midtime.
%   bootAnchor.ub                       : Upper bounds for minimization of size 1 x number of rate constants * number of regions
%   bootAnchor.lb                       : Lower bounds for minimization of size 1 x number of rate constants * number of regions
%   tacVascCorr.vascCorrTacs            : Vascularly corrected time activity curves (TACs) for brain regions in "regions" in the same order as m x r array where m is length(midtime) and r is length(regions)
%   fitVascMeanTac.timeVascFitUniform   : the full scan time values for VascTacFitUniform with uniform sampling.
%   weights                             :vector of n weights for TAC frames (n,1)

out=in;
%For now, STARE has only been run with equal weighting of target TACs. Preliminary validation showed no advantage of weightin by
%region (e.g., there was not preliminary effect of region size when
%considering different regions acting as sources.
out.regionWeights=ones(1,length(in.regions)-1);
%Set up simulannealbnd options
opts = optimoptions('simulannealbnd',...
    'PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplotstopping,@saplottemperature},...
    'InitialTemperature',100,... %matlab default=100 %%%%%%% 'TemperatureFcn','temperatureexp',... %default=temperatureexp %%%%%%% %'AnnealingFcn','annealingfast',... %default=annealingfast &&&&&&&&&     'AcceptanceFcn','acceptancesa',... %default=acceptancesa
    'ReannealInterval',100,... %default=100
    'HybridFcn',[],... %can add fminsearch, patternsearch, fminunc, fmincon as an additional min step (also specify 'HybridInterval','end')
    'FunctionTolerance',1e-6,... %default=1e-6
    'MaxIterations',Inf,... %default=Inf
    'Display','final'); %default=final
out.minimizeCostFunc.opts=opts;
%Get free parameter initial guasses for minimization
out.minimizeCostFunc.x0=in.bootAnchor.lb + (in.bootAnchor.ub - in.bootAnchor.lb).*rand;
%Set up a reduced set of inputs to enter into simulated annealing
saIn={};
saIn.regions=in.regions;
saIn.subOutPath=in.subOutPath;
saIn.vascCorrPerc=in.vascCorrPerc;
saIn.x0=out.minimizeCostFunc.x0;
saIn.midtime=in.midtime;
saIn.tacs=in.tacVascCorr.vascCorrTacs;
saIn.weights=in.weights;
saIn.regionWeights=out.regionWeights;
saIn.lb=in.bootAnchor.lb;
saIn.ub=in.bootAnchor.ub;
saIn.bootKiKsDensPeak=in.bootAnchor.bootKiKsDensPeak;
saIn.timeVascFitUniform=in.fitVascMeanTac.timeVascFitUniform;

singleSourceResults={};
%Initiate parallel pool
parpool(length(in.regions));
tic
%Iterate minimization over each region acting as the source
parfor s=1:length(in.regions)
    %Create a temporary structure to be used within the parfor loop.
    tempIn = saIn;
    tempIn.sourceInd=s;
    tempIn.sourceTac=tempIn.tacs(:,tempIn.sourceInd);
    tempIn.sourceName=tempIn.regions{s};

    tempIn.targetTacs=tempIn.tacs;
    tempIn.targetTacs(:,tempIn.sourceInd)=[]; %Eliminate source TAC, so just target TACs leftover
    tempIn.targetNames=tempIn.regions;
    tempIn.targetNames(tempIn.sourceInd)=[];
    %Must use anonymous function that tells simulannealbnd that x are the parameters optimized.
    objFunc=@(x) costFunc(x,tempIn);
    %Run simulated annealing for this source region
    [finalRateConsts,finalCost,exitFlag,saOut]=simulannealbnd(objFunc,tempIn.x0,tempIn.lb,tempIn.ub,opts);
    %Apply estimated parameters to get TAC fits
    s2tOut=s2tModel(finalRateConsts,tempIn);
    %Gather outputs from simulannealbnd to save
    tempOut=tempIn;
    tempOut.saExitFlag=exitFlag;
    tempOut.saOutput=saOut;
    tempOut.finalCost=finalCost;
    tempOut.singleSourceFinalRateConsts=finalRateConsts;
    
    tempOut.singleSourceFinalKi=s2tOut.ki;
    tempOut.targetTacFits=s2tOut.targetTacFits;
    tempOut.irf=s2tOut.irf;
    %Get final ratio of TAC fitting term cost to Ki penalty term cost
    tacGofTerm=sum( tempIn.regionWeights.*sum( ( tempIn.weights.*(tempIn.targetTacs - tempOut.targetTacFits).^2 ) ) );
    kiPenaltyTerm=sum(abs(tempOut.singleSourceFinalKi - tempIn.bootKiKsDensPeak));  
    tempOut.ratioTacGof2KiPenalty=tacGofTerm/kiPenaltyTerm;
    %Plot the source TAC, target TACs, and target TAC fits for the single
    %source rotation
    tempPlot={};
    tempPlot.midtime=tempIn.midtime;
    tempPlot.sourceTac=tempIn.sourceTac;
    tempPlot.sourceName=tempIn.sourceName;
    tempPlot.targetTacs=tempIn.targetTacs;
    tempPlot.targetNames=tempIn.targetNames;
    tempPlot.targetTacFits=tempOut.targetTacFits;
    tempPlot.sourceInd=tempIn.sourceInd;
    tempPlot.subOutPath=tempIn.subOutPath;
    plotStareTacFits(tempPlot)

    singleSourceResults{s}=tempOut;
end
toc
delete(gcp('nocreate')) %Close parallel pool session.

%Gather rate constants and Ki estimates from each source rotation
%Create headers to write final RCs and Ki to a cell array that can be
%stored out output to csv
headers=['Subject' strcat([repmat({'K1-'},1,length(in.regions)) repmat({'k2-'},1,length(in.regions)) repmat({'k3-'},1,length(in.regions)) repmat({'Ki-'},1,length(in.regions))],repmat(in.regions,1,4))];
out.finalOutput={};
out.finalOutput=headers;
out.finalOutput{2,1}=in.subject;
out.minimizeCostFunc.singelSourceResults=singleSourceResults;
for s=1:length(in.regions)
    out.minimizeCostFunc.finalSingleSourceRateConstsKi(s,:)=[out.minimizeCostFunc.singelSourceResults{s}.singleSourceFinalRateConsts out.minimizeCostFunc.singelSourceResults{s}.singleSourceFinalKi];
end
%Get final RC and Ki results by averaging across source rotations
out.finalRateConstsKi=mean(out.minimizeCostFunc.finalSingleSourceRateConstsKi);
for i=1:length(out.finalRateConstsKi)
   out.finalOutput{2,i+1}=out.finalRateConstsKi(i);
end
%Write final results to csv
writetable(cell2table(out.finalOutput),fullfile(in.subOutPath,[in.subject '_Final_STARE_RateConstants+Ki.csv']), 'WriteVariableNames', 0);

save(fullfile(in.subOutPath,[in.subject '_STARE_out.mat']),'out')


end

