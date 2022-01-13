function [ weightedResidual,tacFit ] = func2tc(x,in )
%func2tc Runs 2TC compartmental models
%   Required inputs:
%   tracer = string radiotracer e.g. 'FDG'
%   midtime = vector of n midtimes for TAC data (n,1)
%   tac = vector of n TAC values (n,1) - vascular correction must alreay be performed on TACs if it is desired
%   weights = vector of n weights for TAC frames (n,1)
%   timeIf = vector of uniformly sampled m blood sampling time-points (m,1)
%   if = vector of uniformly sampled m blood samples, representing input function (m,1)

K1=x(1);
k2=x(2);
k3=x(3);
if strcmp(in.tracer,'FDG')==1 %If tracer=FDG, then k4=0 and runs 2TC-Irrev
    k4=0;
else
    k4=x(4);
end

%Equations pulled directly from Phelps 1979 for 2TC model
alpha1=( k2 + k3 + k4 - sqrt( (k2 + k3 + k4)^2 - 4*k2*k4 ) )/2;
alpha2=( k2 + k3 + k4 + sqrt( (k2 + k3 + k4)^2 - 4*k2*k4 ) )/2;
IRF = K1/(alpha2 - alpha1)*( (k3 + k4 - alpha1)*exp(-alpha1*in.timeIf) + (alpha2 - k3 - k4)*exp(-alpha2*in.timeIf));
%Convolve impulse response function (IRF) with linearly sampled blood data to get the fitted tissue concentration (C_total)
tacFitUniform=(in.timeIf(2)-in.timeIf(1)).*filter(in.if,1,IRF);

%Sample the resulting total concentration in the tissue back to the
%original midtime sampling so can compare to the TAC data
tacFit=interp1(in.timeIf,tacFitUniform,in.midtime,'pchip');

weightedResidual=in.weights.*(in.tac - tacFit);
end

