function [] = plotStareTacFits(in)
%plotStareTacFits Plots the source and target region time activity curves
%(TACs) and the resulting TAC fits from STARE.
%   Requires a structure as input with the following fields:
%   midtime             : String of midtimes for all dynamic frames considered of size (m,1) where m is # of midtimes
%   sourceTac           : Vascularly corrected time activity curve (TAC) for source region in "regions" in m x 1 array where m is length(midtime)
%   sourceName          : String that is the brain region name of the source.
%   targetTacs          : Vascularly corrected time activity curves (TACs) for target regions in "regions" in the same order as m x r array where m is length(midtime) and r is length(regions)-1
%   targetNames         : A cell array of strings with the target brain region names.
%   targetTacFits       : Fits from the source-to-tissue model of the target TACs. Must be array of the same size as targetTacs.
%   sourceInd           : index for the source region
%   subOutPath          : subject-level output directory. Created in stareAnchoring.m

fig1=figure;
colors=get(gca,'ColorOrder');

maxVal=round(max([max(in.sourceTac),max(max(in.targetTacs)),max(max(in.targetTacFits))]),1);
yMax=maxVal+0.1; %in case round func rounded down to nearest 0.1, add 0.1 to ensure axes don't exclude values. Ceil only works with integers

colorsTargets=colors;
colorsTargets(in.sourceInd,:)=[];

subplot(2,1,1)    
hold on
plot(in.midtime,in.sourceTac,'*','Color',colors(in.sourceInd,:),'LineWidth',2,'MarkerSize',10);
for j=1:length(in.targetNames)
    plot(in.midtime,in.targetTacs(:,j),'o','Color',colorsTargets(j,:),'LineWidth',2,'MarkerSize',10);
end
legend([strcat('Source-',strrep(in.sourceName,'_',' '));strcat('Target-',strrep(in.targetNames,'_',' '))'],'Location','eastoutside')
set(gca,'YLim',[0,yMax],'XLim',[0,in.midtime(end)],'FontSize',20)
xlabel('Minutes')
ylabel('microCi')
title('Source and Target Region TACs')

subplot(2,1,2)
hold on
for j=1:length(in.targetNames)
    plot(in.midtime,in.targetTacs(:,j),'o','Color',colorsTargets(j,:),'LineWidth',2,'MarkerSize',10);
end
for j=1:length(in.targetNames)
    plot(in.midtime,in.targetTacFits(:,j),'Color',colorsTargets(j,:),'LineWidth',2);
end    
legend([strrep(in.targetNames,'_',' ')';strcat('Fit-',strrep(in.targetNames,'_',' '))'],'Location','eastoutside')
set(gca,'YLim',[0,yMax],'XLim',[0,in.midtime(end)],'FontSize',20)
xlabel('Minutes')
ylabel('microCi')
title('Target Region TACs and Fits')

savefig(fullfile(in.subOutPath,['STARE_Target_TAC_Fits_with_'  num2str(in.sourceName) '_as_Source.fig']))
close(fig1);
end

