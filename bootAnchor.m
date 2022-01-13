function [out] = bootAnchor (in)
%bootAnchor Uses decreasing sum of 3 exponentials to fit blood data
% Called from single_souce_METHOD.m
%   Required inputs:
%   tracer                                      : radiotracer used. Only validated for FDG currently.
%   midtime                                     : string of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   vascCorrPerc                                : Percentage from 0 to 100 (0 and 5 used in validaions) for contribution to voxel-wise PET signal from vasculature signal (i.e., if =5, 5% of signal is corrupted from vasclature. For vascular correction of TACs.
%   subOutPath                                  : subject-level output directory. Created in stareAnchoring.m
%   subject                                     : string with subject name (e.g., 'FDG001')
%   weights                                     : Generated frame-wise weights (square root of frame duration). Must be same size as midtime.
%   VascClust.figsPath                          : ouput path for figures
%   petPvc.pvcVascSdTac                         : standard deviation across all voxels in the partial volume corrected vascular TAC for each frame 
%   petPvc.pvcVascMeanTac                       : time activity curve (TAC) from the final PVC'ed vascular cluster 
%   fitVascMeanTac.postPeakMidtimes             : string of midtimes for all dynamic frames **after the peak** considered of size (n,1) where n is # of midtimes
%   fitVascMeanTac.postPeakVascTac              : time activity curve (TAC) from the final PVC'ed vascular cluster **after the peak** of size (n,1)
%   fitVascMeanTac.postpeakWeights              : A set of weights for each entry in postPeakMidtimes of size (n,1)
%   fitVascMeanTac.VascTacFitUniform            : the full scan fit of the pvc'ed vascular mean tac with uniform sampling
%   fitVascMeanTac.timeVascFitUniform           : the full scan time values for VascTacFitUniform with uniform sampling.
%   fitVascMeanTac.peakInd                      : Index of petPvc.pvcVascMeanTac peak

out=in;

postPeakMidtimes = in.fitVascMeanTac.postPeakMidtimes;
postPeakVascTac =  in.fitVascMeanTac.postPeakVascTac;
postPeakWeights = in.fitVascMeanTac.postpeakWeights;
postPeakVascSd = in.petPvc.pvcVascSdTac(in.fitVascMeanTac.peakInd:end);
out.bootAnchor.postPeakVascSd = postPeakVascSd;

%Run bootstrapping 1000 times to generate 1000 curves between 1 standard
%deviation less than the mean signal in the PVC vasculature cluster and 1
%standard deviation above the mean signal. 
out.bootAnchor.bootIters=1000;
for iter=1:out.bootAnchor.bootIters
    lowerBound=postPeakVascTac-postPeakVascSd; 
    upperBound=postPeakVascTac+postPeakVascSd;

    %Randomly generate curve between the upper and lower bounds
    bootCurve=lowerBound+(upperBound-lowerBound).*rand(size(lowerBound));
    %Save curve to array with all 1,000 curves.
    out.bootAnchor.bootCurve(:,iter)=bootCurve;
end

%Determine number of 2TCirr params. If/else set up for easy transition to
%other irreversible radiotracers with 2TCirr, or other tracers with
%reversible kinetics.
if strcmp(in.tracer,'FDG')==1 %%%%%%%% NEED TO UPDATE FOR OTHER TRACERS
    out.num2tcParams=3; %2TCIrrev
else
    error('STARE not yet developed or validated for radiotracers outside of FDG')
end
        
out.boot_RCs=zeros(out.bootAnchor.bootIters,length(in.regions)*out.num2tcParams); %3 for 2TCirrev


% Set up fittype and options.
ft = fittype( 'C1*exp(-lambda1*x)+C2*exp(-lambda2*x)+C3*exp(-lambda3*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Algorithm','Levenberg-Marquardt','MaxFunEvals',60000,'MaxIter',40000,'Weights',postPeakWeights,'Display','Off');

%Set up uniform sampling interval for pre- and post-peak
postPeakTimeUniform=round(postPeakMidtimes(1),1):0.1:postPeakMidtimes(end); 
prePeakTimeUniform=0:.1:round(postPeakMidtimes(1),1)-0.1;
%Combine interpolated pre- and post-peak time
timeBootCurveUniform=[prePeakTimeUniform';postPeakTimeUniform'];

%Interpolate pre-peak raw vasc TAC data
prePeakVascTacUniform=interp1(in.midtime(1:in.fitVascMeanTac.peakInd),in.petPvc.pvcVascMeanTac(1:in.fitVascMeanTac.peakInd),prePeakTimeUniform);
prePeakVascTacUniform(logical(sum([prePeakVascTacUniform<0;isnan(prePeakVascTacUniform)])))=0;

%Set up lsqnonlin options for 2TCirr
lsqnonlinOpts=optimset('LargeScale','off','Display','off');

for iter=1:out.bootAnchor.bootIters
    postPeakBootCurve=out.bootAnchor.bootCurve(:,iter);
       
    %Set up fitting so that we fit the boot curve at the given sampling interval
    [xData, yData] = prepareCurveData( postPeakMidtimes, postPeakBootCurve );
   
    success=0;
    fail=0;
    fprintf('Bootstrap Iter %d of %d.\n',iter,out.bootAnchor.bootIters)
    %Run fitting up to 6 times to allow us to get a successful fit
    while success<1 && fail<6
        try
            % Fit model to data.            
            [fitResult,~,fitOut] = fit(xData,yData,ft,opts);
            if fitOut.exitflag>0 %From Matlab: Positive flags indicate convergence, within tolerances.
                success=success+1;

                %Get fit of uniformly sampled boot curve using fitResult
                postPeakFitUniform=fitResult(postPeakTimeUniform);

                %Combine interpolated pre-peak raw blood data & post-peak fit IF
                bootCurveFitUniform=[prePeakVascTacUniform';postPeakFitUniform];

                %Check if data is uniformly sampled. If not throw error. We need uniform
                %sampling for the convolution in Func2TC.m
                check=[];
                for i=2:length(timeBootCurveUniform)
                    check(i-1)=round(timeBootCurveUniform(i)-timeBootCurveUniform(i-1),3); %Interpollating timeIF seems to allow imprecision at >10 decimal places, so we round to then compare the min and max
                end
                if max(check)~=min(check)
                    error('IF passed to Func2TC is not uniformly sampled.')
                end
            else
                fail=fail+1;
            end
        catch
            fail=fail+1;
        end
    end
    %Now that we have a fitted boot curve, use it as input to two tissue
    %irreversible modeling to give us a rnage of rate constants for each
    %bootstrap iteration.
    if success==1 && ~any(postPeakFitUniform<0) %Don't use the bootstrap sample if anywhere after the peak of the fitted curve dips below zero
        out.bootAnchor.bootCurveFits{iter}=[timeBootCurveUniform bootCurveFitUniform];
        %Set 2TCirr upper and lower bounds for K1, k2, k3 to [0,1]
        lb=zeros(1,out.num2tcParams); 
        ub=ones(1,out.num2tcParams);
        %Set up structure to input into Func2TC in lsqnonlin for 2TCirr
        in2tc={};
        in2tc.tracer=in.tracer;
        in2tc.midtime=in.midtime;
        in2tc.timeIf=timeBootCurveUniform;
        in2tc.if=bootCurveFitUniform; %Treating bootstrapped vasc signal curve as input to 2TCirr modeling.
        in2tc.weights=in.weights;
        rateConsts=zeros(1,out.num2tcParams*length(in.regions));
        % Get boot curve sampled at midtimes for vascular correction
        bootCurveAtMidtimes=interp1(timeBootCurveUniform,bootCurveFitUniform,in.midtime,'pchip');
        for r=1:length(in.regions)
            %Vascularly correct TAC
            rawTac=in.tacs(:,r);
            vascCorrTac=(1/(1-in.vascCorrPerc/100))*(rawTac - (in.vascCorrPerc/100)*bootCurveAtMidtimes);
            in2tc.tac=vascCorrTac;
            success=0;
            fail=0;
            %Run 2TCirrr. Give it 6 tries to succeed.
            while success<1 && fail<6
                try
                    x0=lb+(ub-lb).*rand; %Generate random guesses for rate constants within upper and lower bounds
                    [xOpt,~,~,exitFlag]= lsqnonlin('func2tc',x0,lb,ub,lsqnonlinOpts,in2tc);
                    if exitFlag>0 %Could be 1: Function converged to a solution x. 2: Change in x was less than the specified tolerance. 3: Change in the residual was less than the specified tolerance. 4: Magnitude of search direction was smaller than the specified tolerance.
                        success=success+1;
                        %K1
                        rateConsts(r)=real(xOpt(1));
                        %k2
                        rateConsts(r+length(in.regions))=real(xOpt(2));
                        %k3
                        rateConsts(r+length(in.regions)*2)=real(xOpt(3));
                    else
                        fail=fail+1;
                    end
                catch
                    fail=fail+1;
                end
            end
        end
            if any(round(rateConsts,2)==0) || any(round(rateConsts,1)==1) %Throw away runs where 2TC irrev yields 0 or 1 values
                rateConsts(:)=nan;
            end
    else
        %If we can't get a fitResult for the boot curve at hand, set rate constants to NaN and try next iteration
        rateConsts(:)=nan;
        out.bootAnchor.bootCurveFits{iter}=[nan nan];
    end
out.bootAnchor.bootRateConsts(iter,:)=rateConsts;
end

% % % %Now that we have the 1,000 sets of bootstrapped rate constants, use them
% % % %to generate probability density functions, then take the full width at
% % % %half max of that density function to get the range of free parameters in
% % % %STARE and the penalty in the cost function


%%%Assign STARE upper and lower bounds as either side of ksdensity FWHM of
%%%%bootstrap samples
%BS Bounds for K1, k2, and k3 (for constraining STARE search spaceout.bootAnchor.lb=zeros(1,out.num2tcParams); 
out.bootAnchor.ub=ones(1,out.num2tcParams);
data.bootstrap_ksdens_peak=out.bootAnchor.ub;
for k=1:out.num2tcParams*length(in.regions)
    [f,xI]=ksdensity(out.bootAnchor.bootRateConsts(:,k));
    xINew=linspace(min(out.bootAnchor.bootRateConsts(:,k)),max(out.bootAnchor.bootRateConsts(:,k)),1000); %xi & f are size(bootstat,1) -- want more frequent + constant sampling independent of # of bootstrap samples
    fNew=interp1(xI,f,xINew);
    out.bootAnchor.lb(k)=xINew(find(fNew>=max(fNew)/2,1,'first')); %Finds rate constant value at left (lower) side of the ksdensity FWHM
    out.bootAnchor.ub(k)=xINew(find(fNew>=max(fNew)/2,1,'last')); %Finds rate constant value at right (upper) side of the ksdensity FWHM
    [~,i]=max(fNew);
    out.bootAnchor.bootKsDensPeak(k)=xINew(i);
end
out.bootAnchor.stareBounds=[out.bootAnchor.ub;out.bootAnchor.lb];
%
%BS KS density peak for Ki (for penalty term in cost function)
bootK1=out.bootAnchor.bootRateConsts(:,1:length(in.regions));
bootK2=out.bootAnchor.bootRateConsts(:,length(in.regions)+1:length(in.regions)*2);
bootK3=out.bootAnchor.bootRateConsts(:,length(in.regions)*2+1:length(in.regions)*3);
out.bootAnchor.bootKis=bootK1.*bootK3./(bootK2+bootK3);
%BS Bounds for Ki
for k=1:length(in.regions)
    [f,xI]=ksdensity(out.bootAnchor.bootKis(:,k));
    xINew=linspace(min(out.bootAnchor.bootKis(:,k)),max(out.bootAnchor.bootKis(:,k)),1000); %xi & f are size(bootstat,1) -- want more frequent + constant sampling independent of # of bootstrap samples
    fNew=interp1(xI,f,xINew);
    [~,i]=max(fNew);
    out.bootAnchor.bootKiKsDensPeak(k)=xINew(i);
end

plotIn={};
plotIn.regions=in.regions;
plotIn.bootRateConsts=out.bootAnchor.bootRateConsts;
plotIn.bootKis=out.bootAnchor.bootKis;
plotIn.subject=in.subject;
plotIn.figsPath=in.vascClust.figsPath;
plotIn.bounds=out.bootAnchor.stareBounds;
plotIn.bootCurve=out.bootAnchor.bootCurve;
plotIn.postPeakMidtimes=postPeakMidtimes;
plotIn.bootCurveFits=out.bootAnchor.bootCurveFits;
plotIn.midtime=in.midtime;
plotIn.pvcVascMeanTac=in.petPvc.pvcVascMeanTac;
plotIn.pvcVascTacFit=in.fitVascMeanTac.VascTacFitUniform;
plotIn.timeVascFitUniform=in.fitVascMeanTac.timeVascFitUniform;
plotBootAnchor(plotIn);

save(fullfile(in.subOutPath,[in.subject '_STARE_out.mat']),'out')

end

