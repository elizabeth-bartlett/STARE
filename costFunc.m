function [ cost ] = costFunc( x,in )
%costFunc Executes STARE objective function
% See Bartlett et al. 2022 in NeuroImage for cost function.
% Reguired Inputs:
%   regionWeights
%   weights                             : Generated frame-wise weights (square root of frame duration). Must be same size as midtime.
%   targetTacs                          : Vascularly corrected time activity curves (TACs) for target regions in "regions" in the same order as m x r array where m is length(midtime) and r is length(regions)-1
%   sourceTac                           : Vascularly corrected time activity curve (TAC) for source region in "regions" in m x 1 array where m is length(midtime)
%   sourceInd                           : index for the source region
%   regions                             : Brain region names as strings in a cell array that will be quantified in STARE.
%   midtime                             : string of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   bootKiKsDensPeak                    : The peak of the probability density estimate for each region across bootstrapping intervals for Ki of size 1 x the number of regions
%   timeVascFitUniform                  : the full scan time values for VascTacFitUniform with uniform sampling.

%   Fit target TACs relative to the given source TAC.
[s2tOut]=s2tModel(x,in);

%Impose penalty for Ki
kiPenalty=sum(abs(s2tOut.ki - in.bootKiKsDensPeak));  

%Compute cost (y) for objective function that gives squared difference between target TACs and TAC fits across all target TACs fit as a function of the common source with weights for regions and frames and the Ki penalty term.
cost=sum( in.regionWeights.*sum( ( in.weights.*(in.targetTacs-s2tOut.targetTacFits).^2 ) ) ) + kiPenalty;
end