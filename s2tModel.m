function [ out ] = s2tModel (x, in )
%s2tModel Fits the TACs according to the source-to-target tissue model
% See Bartlett et al. 2022 in NeuroImage for full s2t derivation in the
% appendix.
% Reguired Inputs:
%   x                   : Vector of free parameter estimates
%   sourceInd           : index for the source region
%   regions             : Brain region names as strings in a cell array that will be quantified in STARE.
%   midtime             : string of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   timeVascFitUniform  : the full scan time values for VascTacFitUniform with uniform sampling.
%   sourceTac           : the time activity curve (TAC for the source region of size m x 1 where m is length(midtime)

out=in;

%SOURCE TAC ---- Extract parameter estimates for the source region
    k1S=x(in.sourceInd);
    k2Se=x(in.sourceInd+length(in.regions));
    k3S=x(in.sourceInd+length(in.regions)*2);
%TARGET TAC ---- Extract parameter estimates for the target regions. 2nd line of each rate constant eliminates source region parameter from target parameter list
    %First get K1, k2, and k3 for all regions (including source)
    k1T=x(1:length(in.regions));
    k2T=x(length(in.regions)+1:length(in.regions)*2);
    k3T=x(length(in.regions)*2+1:end);    
    %Before eliminate source rate constants (to only have target rate
    %constants as we need for the rest of the script, compute Ki for both
    %the source and targets at once.
    out.ki=(k1T.*k3T)./(k2T+k3T);
    %Now remove source rate constants to only leave target rate constants
    k1T(in.sourceInd)=[];
    k2T(in.sourceInd)=[];
    k3T(in.sourceInd)=[];

%%%%% Expressions used to express the target regions in terms of the source
%%%%% region. Refer to STARE Appendix for these equations and fariable
%%%%% naming

%Single valued p-,q-,rS (for the source). Use repmat to repeat the single value
%to account for each of the target regions
pS=repmat(k2Se/(k2Se+k3S),1,length(in.regions)-1);
qS=repmat(k3S/(k2Se+k3S),1,length(in.regions)-1);
rS=repmat(k2Se+k3S,1,length(in.regions)-1);

pT = k2T./(k2T+k3T);
qT = k3T./(k2T+k3T);
rT = k2T+k3T;

alpha = (qT.*rT)./(pT+qT) + rS - (qS.*rS)./(pS+qS) - rT;
beta = (qT.*rT.*rS)./(pT+qT) - (qS.*rT.*rS)./(pS+qS);
gamma = (qS.*rS)./(pS+qS) + rT;
omega = (qS.*rT.*rS)./(pS+qS);

nu = ( -gamma + sqrt(gamma.^2 - 4*omega) )/2;
epsilon = ( -gamma - sqrt(gamma.^2 - 4*omega) )/2;

L = alpha - (beta + alpha.*epsilon)./(epsilon-nu);
M = (beta + alpha.*epsilon)./(epsilon-nu);

%Must do fitting with uniformly sampled data
sourceTacUniform = interp1(in.midtime,in.sourceTac,in.timeVascFitUniform,'pchip');
%CHeck to ensure sampling is uniform
check=[];
for i=2:length(in.timeVascFitUniform)
    check(i-1)=round(in.timeVascFitUniform(i)-in.timeVascFitUniform(i-1),3); %Interpollating time seems to allow imprecision at >10 decimal places, so we round to then compare the min and max
end
if max(check)~=min(check)
    error('IF passed to Func2TC is not uniformly sampled.')
end
samplingInterval=max(check);

for target=1:length(in.regions)-1 %Total regions - 1 = # of target regions being related to the source
   %Again, reference STARE appendix for equations.
   zeta = L(target)*exp(nu(target)*in.timeVascFitUniform) + M(target)*exp(epsilon(target)*in.timeVascFitUniform);
   %Compute target TAC fit as a function of the source TAC. Source-to-target tissue model
   targetTacFitUniform = (k1T(target)/k1S)*samplingInterval*filter(sourceTacUniform,1,zeta) + (k1T(target)/k1S)*sourceTacUniform;
   
   out.irf(:,target)=interp1(in.timeVascFitUniform,zeta,in.midtime,'pchip');
   out.targetTacFits(:,target)=interp1(in.timeVascFitUniform,targetTacFitUniform,in.midtime,'pchip');
end

end
