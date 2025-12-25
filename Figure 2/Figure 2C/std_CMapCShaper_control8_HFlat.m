clear all;clc;close all;

load(['..\Figure 2B\CellVariationTime_CMapCShaper_control8.mat']);CellVariCMap=CellVari;
load(['..\Figure 2B\CellVariationTime_CShaperCMap_control8.mat']);CellVariCShaper=CellVari;
TimeList=-10:5:200;

CMapColor=[32,118,203]./255;CShaperColor=[99,185,33]./255;

figure
ShapeDescriptorName='Hayakawa Flatness Ratio';
% CMap
[xTime,yDMean]=Boxplot(CellVariCMap,ShapeDescriptorName,TimeList);
plot(xTime,yDMean,'.','MarkerSize',10,'Color',CMapColor);hold on;
h1=plot(xTime,yDMean,'-','LineWidth',1.5,'Color',CMapColor);hold on;

% CShaper
[xTime,yDMean]=Boxplot(CellVariCShaper,ShapeDescriptorName,TimeList);
plot(xTime,yDMean,'.','MarkerSize',10,'Color',CShaperColor);hold on;
h2=plot(xTime,yDMean,'-','LineWidth',1.5,'Color',CShaperColor);hold on;

plot([150,150],[0.08,1],'--','LineWidth',1.5,'Color',[169,169,169]./255);
legend([h1,h2],'under natural condition','with mechanical compression','Location','southeast','FontSize',13,'FontName','Arial');
legend('boxoff')

axis square;
ylabel({'Variation of ',ShapeDescriptorName});
yticks(0.05:0.05:0.3);
xticks(0:50:250);xlabel('Developmental time (min)');
xlim([-15,205]);ylim([0.05,0.21]);
set(gca,'FontSize',15,'FontName','Arial');
set(gca,'unit','centimeters','position',[1.5,1.6,12,8]);
set(gcf,'unit','centimeters','position',[10,5,14,10]);

function [xTime,yDMean]=Boxplot(CellVari,ShapeDescriptorName,TimeList)
Row=find(strcmp(CellVari(1,:),ShapeDescriptorName));
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