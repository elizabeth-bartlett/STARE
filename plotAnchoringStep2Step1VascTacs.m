function [] = plotAnchoringStep2Step1VascTacs( in )
%plotAnchoringStep2Step1VascTacs Creates matlab .fig with both the step 1 
%and step 2 vascular clusters. Function nested in vascClustering.
%   Requires a structure as input with the following fields:
%   step1VascMeanTac        : time activity curve (TAC) from the final step 1 vascular cluster
%   step2VascMeanTac        : time activity curve (TAC) from the final step 2 vascular cluster
%   midtime                 : string of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   figsPath                : ouput path
%   subject                 : subject ID

figure
subplot(1,3,1)
plot(in.midtime,in.step1VascMeanTac,'LineWidth',3)
hold on
plot(in.midtime,in.step2VascMeanTac,'LineWidth',3)
xlabel('Minutes')
ylabel('MicroCuries')
title([in.subject ' - Whole Scan'])
set(gca,'YLim',[0 inf],'XLim',[0 in.midtime(end)],'FontSize',15)

subplot(1,3,2)
plot(in.midtime,in.step1VascMeanTac,'LineWidth',3)
hold on
plot(in.midtime,in.step2VascMeanTac,'LineWidth',3)
legend('Step 1 Mean Vasc TAC','Step 2 Mean Vasc TAC')
xlabel('Minutes')
ylabel('MicroCuries')
title([in.subject ' - Early-Scan'])
set(gca,'YLim',[0 inf],'XLim',[0 5],'FontSize',15)

subplot(1,3,3)
plot(in.midtime,in.step1VascMeanTac,'LineWidth',3)
hold on
plot(in.midtime,in.step2VascMeanTac,'LineWidth',3)
xlabel('Minutes')
ylabel('MicroCuries')
title([in.subject ' - Late-Scan'])
set(gca,'XLim',[10 in.midtime(end)],'FontSize',15)

savefig(fullfile(in.figsPath,'Step2+Step1_vasc-TACs.fig'))
close all
end

