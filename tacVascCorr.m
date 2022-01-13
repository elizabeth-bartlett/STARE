function [ out ] = tacVascCorr( in )
%tacVascCorr Vascularly corrects TACs based on the fitted mean signal from the partial
%volume corrected vasculature cluster
%   Required inputs:
%   regions                         : Brain region names as strings in a cell array that will be quantified in STARE.
%   tacs                            : Time activity curves (TACs) for brain regions in "regions" in the same order as m x r array where m is length(midtime) and r is length(regions)
%   vascCorrPerc                    : Percentage from 0 to 100 (0 and 5 used in validaions) for contribution to voxel-wise PET signal from vasculature signal (i.e., if =5, 5% of signal is corrupted from vasclature. For vascular correction of TACs.
%   fitVascMeanTac.pvcVascTacFit    : Fitted mean signal from PVCed vasculature cluster sampled at midtimes
%   subOutPath                      : subject-level output directory. Created in stareAnchoring.m
%   subject                         : string with subject name (e.g., 'FDG001')

out=in;

out.tacVascCorr.vascCorrTacs=zeros(size(in.tacs));
for r=1:length(in.regions)
    out.tacVascCorr.vascCorrTacs(:,r)=(1/(1-in.vascCorrPerc/100))*(in.tacs(:,r) - (in.vascCorrPerc/100)*in.fitVascMeanTac.pvcVascTacFit);
end

save(fullfile(in.subOutPath,[in.subject '_STARE_out.mat']),'out')
