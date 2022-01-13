function [] = plotAnchoringStep1VascTacs( in )
%plotAnchoringStep1VascTacs Creates matlab .fig with each of the identified
%vascular TACs from each k-means iteration (varying k). Function nested in
%vascClustering.
%   Requires a structure as input with the following fields
%   step1MeanTacs           : time activity curves (TACs) from each cluster for each iteration of k-means
%   midtime                 : string of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%	ks                      : Vector of k clusters run in each k-means iteration.
%	vascClustInds           : Indices, referring to step1MeanTacs, for the vascular TAC from each k-means run
%	step1VascMeanTac        : The final Step 1 vascular mean TAC
%   step1OptimalNumClus     : The optimal cluster number (k-means run) that produced the final vascular TAC
%   figsPath                : ouput path

figure
hold on
for k=1:length(in.ks)
    if ~isnan(in.vascClustInds(k))
        plot(in.midtime,in.step1MeanTacs(:,in.vascClustInds(k),k),'.-','MarkerSize',30,'LineWidth',2)
    end
end
plot(in.midtime,in.step1VascMeanTac,'k.-','MarkerSize',30,'LineWidth',3)
legend([strcat('k=',string(in.ks(~isnan(in.vascClustInds)))) strcat('Selected Vasc TAC (k=',string(in.step1OptimalNumClus),')')],'Location','eastoutside')
xlabel('Minutes')
ylabel('Activity/cc')
set(gca,'FontSize',20)
title(strcat('Optimal Vascular TACs:',num2str(in.ks(1)),'-',num2str(in.ks(end)),' k-means clusters'))

savefig(fullfile(in.figsPath,'Step1_potential-vasc-clusts.fig'))
% close(fig);
end

