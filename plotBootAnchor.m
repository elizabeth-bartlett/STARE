function [] = plotBootAnchor(in)
%plotBootAnchor Plots the bootstrapped curves and resulting histograms used
%for STARE anchoring
%   Requires a structure as input with the following fields:
%   regions                             : Brain region names as strings in a cell array that will be quantified in STARE.
%   subject                             : string with subject name (e.g., 'FDG001')
%   midtime                             : string of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   VascClust.figsPath                  : ouput path
%   fitVascMeanTac.postPeakMidtimes     : string of midtimes for all dynamic frames **after the peak** considered of size (n,1) where n is # of midtimes
%   petPvc.pvcVascMeanTac               : time activity curve (TAC) from the final PVC'ed vascular cluster
%   fitVascMeanTac.VascTacFitUniform    : the full scan fit of the pvc'ed vascular mean tac with uniform sampling
%   fitVascMeanTac.timeVascFitUniform   : the full scan time values for VascTacFitUniform with uniform sampling.
%   bootAnchor.bootRateConsts           : A b x r array of rate constants estimated from passing bootstrapped vasculature curve through 2TCirr, where b is the number of bootstrap iterations and r is the number of regions * the number of rate constants (i.e., 3 for 2TCirr with FDG).
%   bootAnchor.bootKis                  : A b x k array of Ki estimates from passing bootstrapped vasculature curve through 2TCirr, where b is the number of bootstrap iterations and r is the number of regions.
%   bootAnchor.bounds                   : 2 x r array of rate constant bounds (all free parameters) where ub is row 1 and lb is row 2 and r is the number of regions * the number of rate constants (i.e., 3 for 2TCirr with FDG).
%   bootAnchor.bootCurve                : A m x b array of bootstrapped curves where m is length(fitVascMeanTac.postPeakMidtimes) and b is the number of bootstrap iterations.
%   bootAnchor.bootCurveFits            : A cell array of length b (the number of bootstrap iterations, where each iteration field (i.e., booCurveFits{1}) contains a s x 2 array where s is length(fitVascMeanTac.timeVascFitUniform), the 1st column is the timee values and the 2nd value is the activity values of the fitted bootstrap curve.

f1=figure;
f2=figure;
f3=figure;
f4=figure;

for i=1:length(in.regions)*4
    if i<=length(in.regions)
        set(0,'CurrentFigure',f1,'DefaultLineLineWidth',3)
        subplot(2,3,i)
    elseif i>length(in.regions) && i<=length(in.regions)*2
        set(0,'CurrentFigure',f2,'DefaultLineLineWidth',3)
        subplot(2,3,i-length(in.regions))
    elseif i>length(in.regions)*2 && i<=length(in.regions)*3
        set(0,'CurrentFigure',f3,'DefaultLineLineWidth',3)
        subplot(2,3,i-length(in.regions)*2)
    elseif i>length(in.regions)*3 && i<=length(in.regions)*4
        set(0,'CurrentFigure',f4,'DefaultLineLineWidth',3)
        subplot(2,3,i-length(in.regions)*3)
    end
    
    if i<=length(in.regions)*3
        [f,xi]=ksdensity(in.bootRateConsts(:,i));
        hold on
        h=hist(in.bootRateConsts(:,i),50);
        hist(in.bootRateConsts(:,i),50)
    elseif i>length(in.regions)*3
        [f,xi]=ksdensity(in.bootKis(:,i-length(in.regions)*3));
        hold on
        h=hist(in.bootKis(:,i-length(in.regions)*3),50);
        hist(in.bootKis(:,i-length(in.regions)*3),50)
    end
    
    ylabel('# of bootstrap iters')
    if i<=length(in.regions)
        title(sprintf('%s K_1 %s',in.subject,in.regions{i}))
        xlabel('Micro-parameter estimate')
    elseif i>length(in.regions) && i<=length(in.regions)*2
        title(sprintf('%s k_2 %s',in.subject,in.regions{i-length(in.regions)}))
        xlabel('Micro-parameter estimate')
    elseif i>length(in.regions)*2 && i<=length(in.regions)*3
        title(sprintf('%s k_3 %s',in.subject,in.regions{i-length(in.regions)*2}))
        xlabel('Micro-parameter estimate')
    elseif i>length(in.regions)*3 && i<=length(in.regions)*4
        title(sprintf('%s K_i %s',in.subject,in.regions{i-length(in.regions)*3}))
        xlabel('Macro-parameter estimate')
    end
    set(gca,'FontSize',15)
    
    n=max(h)*(f-min(f))./(max(f)-min(f));
    findpeaks(n,xi)
    set(gca,'YLim',[0 Inf],'XLim',[-Inf Inf])
    
    if i<=length(in.regions)*3
        plot([in.bounds(1,i),in.bounds(2,i)],[max(h)/2 max(h)/2],'c')
        plot([in.bounds(1,i),in.bounds(1,i)],[0 max(h)/2],'c')
        plot([in.bounds(2,i),in.bounds(2,i)],[0 max(h)/2],'c')
        legend('bs"ed hist','ksdens','peak','bs"ed bounds')
    else
        legend('bs"ed hist','ksdens','peak')
    end
    
end

%Boot Vasc Curve plot
f5=figure;
title(sprintf('%s Curves from Bootstrapped Vasculature Cluster',in.subject))
    xlabel('Minutes')
    ylabel('Activity')
    set(gca,'FontSize',20)
    hold on
for i=1:length(in.bootCurveFits)  
    plot(in.bootCurveFits{i}(:,1),in.bootCurveFits{i}(:,2),'LineWidth',2)
    plot(in.postPeakMidtimes,in.bootCurve(:,i),'o','MarkerSize',8,'LineWidth',1)
end
    plot(in.midtime,in.pvcVascMeanTac,'*','MarkerSize',10,'LineWidth',2)
    plot(in.timeVascFitUniform,in.pvcVascTacFit,'--k','LineWidth',3)

savefig(f1,fullfile(in.figsPath,'boot_K1.fig'));
savefig(f2,fullfile(in.figsPath,'boot_k2.fig'));
savefig(f3,fullfile(in.figsPath,'boot_k3.fig'));
savefig(f4,fullfile(in.figsPath,'boot_Ki.fig'));
savefig(f5,fullfile(in.figsPath,'boot_curves.fig'));

close all
end
