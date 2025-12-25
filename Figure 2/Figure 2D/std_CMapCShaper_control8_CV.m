clear all;clc;close all;

load(['.\CellVariationTime_CMapCShaper_ContactCV_control8.mat']);CellVariWT=CellVari;
load(['.\CellVariationTime_CShaperCMap_ContactCV_control8.mat']);CellVariCShaper=CellVari;
% time interval
TimeList=0:5:240;

CMapColor=[32,118,203]./255;CShaperColor=[99,185,33]./255;

figure
CVName='CV';
% WT
[xTime,yDMeanWT]=Boxplot(CellVariWT,CVName,TimeList);
plot(xTime,yDMeanWT,'.','MarkerSize',10,'Color',CMapColor);hold on;
h1=plot(xTime,yDMeanWT,'-','LineWidth',1.5,'Color',CMapColor);hold on;

% lag-1
[xTime,yDMeanCShaper]=Boxplot(CellVariCShaper,CVName,TimeList);
plot(xTime,yDMeanCShaper,'.','MarkerSize',10,'Color',CShaperColor);hold on;
h2=plot(xTime,yDMeanCShaper,'-','LineWidth',1.5,'Color',CShaperColor);hold on;

plot([150,150],[0.4,10],'--','LineWidth',1.5,'Color',[169,169,169]./255);
legend([h1,h2],'under natural condition','with mechanical compression','Location','southeast','FontSize',13,'FontName','Arial');
legend('boxoff');

axis square;
ylabel({'Variation of','cell-cell contact area'});

xticks(0:50:250);xlabel('Developmental time (min)');
yticks(0:0.4:5);
xlim([-15,210]);ylim([-0.1,2]);
set(gca,'FontSize',15,'FontName','Arial');
set(gca,'unit','centimeters','position',[1.5,1.6,12,8]);
set(gcf,'unit','centimeters','position',[10,5,14,10]);

function [xTime,yDMean]=Boxplot(CellVari,PVName,TimeList)
Row=find(strcmp(CellVari(1,:),PVName));
TimeData=cell2mat(CellVari(2:end,2));
DescriptorData=cell2mat(CellVari(2:end,Row));
% classification
xTime=[];yDMean=[];yDStd=[];
for Timebox=1:length(TimeList)-1
    Row=find(TimeData>=TimeList(Timebox) & TimeData<TimeList(Timebox+1));
    TimeTemp=TimeData(Row);DescriptorTemp=DescriptorData(Row);
    xTime=[xTime;(TimeList(Timebox)+TimeList(Timebox))./2];
    yDMean=[yDMean;mean(DescriptorTemp)];
end
end
