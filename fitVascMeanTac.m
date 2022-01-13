function [out] = fitVascMeanTac (in)
%fitVascMeanTac Uses decreasing sum of 3 exponentials to fit mean time activity
%cufrve from vasculature, with interpollation from time zero to the peak.
%   Required inputs:
%   midtime                 : vector of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   petPvc.pvcVascMeanTac   : time activity curve (TAC) from the final PVC'ed vascular cluster
%   subOutPath              : subject-level output directory. Created in stareAnchoring.m
%   subject                 : string with subject name (e.g., 'FDG001')

out=in;

%Get weights for fitting (square root of frame duration)
endTimeFrame=filter(2,[1 1],in.midtime);% Generate the end times from the frame midtimes 
frameDurs=endTimeFrame(1:end,1)-[0;endTimeFrame(1:end-1,1)];
out.weights=real(sqrt(frameDurs));

%Get midtime and vasc TAC data post-peak   
[~,peakInd]=max(in.petPvc.pvcVascMeanTac);
out.fitVascMeanTac.peakInd=peakInd;
out.fitVascMeanTac.postPeakMidtimes=in.midtime(peakInd:end);
out.fitVascMeanTac.postPeakVascTac=in.petPvc.pvcVascMeanTac(peakInd:end);
out.fitVascMeanTac.postpeakWeights=out.weights(peakInd:end);

%Set up fitting
[xData, yData] = prepareCurveData(out.fitVascMeanTac.postPeakMidtimes,out.fitVascMeanTac.postPeakVascTac);
%Set up fittype and options.
ft = fittype( 'C1*exp(-lambda1*x)+C2*exp(-lambda2*x)+C3*exp(-lambda3*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Algorithm','Levenberg-Marquardt','MaxFunEvals',60000,'MaxIter',40000,'Weights',out.fitVascMeanTac.postpeakWeights,'Display','Off');

success=0;
fail=0;

sqWeightedResiduals=[];
fitResults={};
%Run fitting up to 10000 times to allow us to get 10 successful fitting runs
while success<10 && fail<10000
    try
        %Fit model to data.
        [fitResult,~,fitOut] = fit(xData,yData,ft,opts);
        warning('off','last');
        if fitOut.exitflag>=0 %From Matlab: Positive flags indicate convergence, within tolerances. 0 indicates hit max # of evals.
            success=success+1;
            fprintf('Vasculature TAC Fitting Success %d.\n',success)
            fitResults{success}=fitResult;
            postPeakVascTacFit=fitResult(out.fitVascMeanTac.postPeakMidtimes); %Apply fit parameters estimated in fitresult to the post-peak midtimes
            %Calculate weighted error
            sqWeightedResiduals(success)=sum( ( out.fitVascMeanTac.postpeakWeights.*(yData - postPeakVascTacFit).^2 ) );
        else
            fail=fail+1;
        end
    catch
        fail=fail+1;
    end
end
out.fitVascMeanTac.fitResults=fitResults;
out.fitVascMeanTac.sqWeightedResiduals=sqWeightedResiduals;
if success==10
    fprintf('-------------------------------\nOptimal Vasculature TAC fit using found considering the best fit out of 10 optimizations.\n-------------------------------\n')

    %Set up uniform sampling interval for pre- and post-peak
    postPeakTimeUniform=round(out.fitVascMeanTac.postPeakMidtimes(1),1):0.1:out.fitVascMeanTac.postPeakMidtimes(end); 
    prePeakTimeUniform=0:.1:round(postPeakTimeUniform(1),1)-0.1;

    %Find the fit with the smallest weighted sum os squared errors
    [~,out.fitVascMeanTac.indBestFit] = min(sqWeightedResiduals);
    %Extract the optimal fitResult from all fitResults
    vascTacFitResult=fitResults{out.fitVascMeanTac.indBestFit};
    out.fitVascMeanTac.finalFitResult=vascTacFitResult;
    %Get  uniformly sampled fit using optimal fitResult
    postPeakVascTacFitUniform=vascTacFitResult(postPeakTimeUniform);

    %Interpolate pre-peak raw blood data.
    prePeakVascTacUniform=interp1(in.midtime(1:peakInd),in.petPvc.pvcVascMeanTac(1:peakInd),prePeakTimeUniform);
    prePeakVascTacUniform(logical(sum([prePeakVascTacUniform<0;isnan(prePeakVascTacUniform)])))=0;

    %Combine interpolated pre-peak raw vasculature TAC with post-peak fit
    %to get final uniformly sampled fit
    out.fitVascMeanTac.timeVascFitUniform=[prePeakTimeUniform';postPeakTimeUniform'];
    out.fitVascMeanTac.VascTacFitUniform=[prePeakVascTacUniform';postPeakVascTacFitUniform];
    
    %Lastly, resample back to the PET midtimes for later vascular correction of
    %the TACs
    out.fitVascMeanTac.pvcVascTacFit=interp1(out.fitVascMeanTac.timeVascFitUniform,out.fitVascMeanTac.VascTacFitUniform,in.midtime,'pchip');

elseif success<10
    error('lsqnonlin could not estimate 10 Vasculature TAC fits in 10000 iterations.')
end

save(fullfile(in.subOutPath,[in.subject '_STARE_out.mat']),'out')


end