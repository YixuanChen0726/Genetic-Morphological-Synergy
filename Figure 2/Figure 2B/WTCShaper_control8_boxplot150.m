% comparison of 12 shape descriptor at 15 min
clear all;clc;close all;
TimeList=-10:5:240;

% WT std
load(['.\CellVariationTime_CMapCShaper_control8.mat']);CellVariCMap=CellVari;
% WT:time/cell stage
load(['.\MeanCellNum_WT.mat']);
CMapStage=MeanNumTime;
CMapRelation=cell2mat(CMapStage(2:end,[1,2]));
CellVariCMap{1,15}='Cell Number';
for CellIndex=2:size(CellVariCMap,1)
    Time=CellVariCMap{CellIndex,2};
    [~,n]=min(abs(CMapRelation(:,1)-Time));
    CellVariCMap{CellIndex,15}=CMapRelation(n,2);
end

% CShaper std
load(['.\CellVariationTime_CShaperCMap_control8.mat']);CellVariCShaper=CellVari;
% CShaper:time/cell stage
load(['.\MeanCellNum_CShaper.mat']);
CShaperStage=cell(size(MeanNumTime,1),10);
CShaperStage(:,3:end)=MeanNumTime(:,12:19);
CShaperStage(1,1:2)={'Time List','Mean Cell Number'};CShaperStage(:,1)=MeanNumTime(:,1);
for CellIndex=2:size(MeanNumTime,1)
    Time=[];
    for Column=12:19
        Temp=MeanNumTime{CellIndex,Column};
        if ~isempty(Temp)
            Time=[Time;Temp];
        end
    end
    CShaperStage{CellIndex,2}=mean(Time);
end
CShaperRelation=cell2mat(CShaperStage(2:end,[1,2]));
CellVariCShaper{1,15}='Cell Number';
for CellIndex=2:size(CellVariCShaper,1)
    Time=CellVariCShaper{CellIndex,2};
    [~,n]=min(abs(CShaperRelation(:,1)-Time));
    CellVariCShaper{CellIndex,15}=CShaperRelation(n,2);
end

figure
% time=149.5-150.5
Row=find(cell2mat(CellVariCMap(2:end,2))>=149.5 & cell2mat(CellVariCMap(2:end,2))<=150.5);
VariCMap150={};VariCMap150(1,:)=CellVariCMap(1,:);
VariCMap150=[VariCMap150;CellVariCMap(Row+1,:)];
CMapmean=mean(cell2mat(VariCMap150(2:end,3:14)));
CMapstd=std(cell2mat(VariCMap150(2:end,3:14)));

Row=find(cell2mat(CellVariCShaper(2:end,2))>=149.5 & cell2mat(CellVariCShaper(2:end,2))<=150.5);
VariCShaper150={};VariCShaper150(1,:)=CellVariCShaper(1,:);
VariCShaper150=[VariCShaper150;CellVariCShaper(Row+1,:)];
CShapermean=mean(cell2mat(VariCShaper150(2:end,3:14)));
CShaperstd=std(cell2mat(VariCShaper150(2:end,3:14)));

CMapColor=[32,118,203]./255;CShaperColor=[99,185,33]./255;W=0.8;
% plot
for DesNum=1:12
    patch([DesNum,DesNum+0.4,DesNum+0.4,DesNum],[0,0,CMapmean(DesNum),CMapmean(DesNum)],CMapColor,'EdgeColor','k','LineWidth',W,'FaceAlpha',0.4);hold on;
    patch([DesNum+0.4,DesNum+0.8,DesNum+0.8,DesNum+0.4],[0,0,CShapermean(DesNum),CShapermean(DesNum)],CShaperColor,'EdgeColor','k','LineWidth',W,'FaceAlpha',0.4);hold on;
    plot([DesNum+0.2,DesNum+0.2],[CMapmean(DesNum)-CMapstd(DesNum)./2,CMapmean(DesNum)+CMapstd(DesNum)./2],'Color','k','LineWidth',W);hold on;
    plot([DesNum+0.13,DesNum+0.27],[CMapmean(DesNum)-CMapstd(DesNum)./2,CMapmean(DesNum)-CMapstd(DesNum)./2],'Color','k','LineWidth',W);hold on;
    plot([DesNum+0.13,DesNum+0.27],[CMapmean(DesNum)+CMapstd(DesNum)./2,CMapmean(DesNum)+CMapstd(DesNum)./2],'Color','k','LineWidth',W);hold on;
    plot([DesNum+0.6,DesNum+0.6],[CShapermean(DesNum)-CShaperstd(DesNum)./2,CShapermean(DesNum)+CShaperstd(DesNum)./2],'Color','k','LineWidth',W);hold on;
    plot([DesNum+0.53,DesNum+0.67],[CShapermean(DesNum)-CShaperstd(DesNum)./2,CShapermean(DesNum)-CShaperstd(DesNum)./2],'Color','k','LineWidth',W);hold on;
    plot([DesNum+0.53,DesNum+0.67],[CShapermean(DesNum)+CShaperstd(DesNum)./2,CShapermean(DesNum)+CShaperstd(DesNum)./2],'Color','k','LineWidth',W);hold on;
end
x1=1.3;x2=1.8;
y1=0.23;y2=0.19;dy=0.02;
patch([x1,x2,x2,x1],[y1-dy y1-dy y1 y1],CMapColor,'EdgeColor','k','LineWidth',W,'FaceAlpha',0.4);hold on;
patch([x1,x2,x2,x1],[y2-dy y2-dy y2 y2],CShaperColor,'EdgeColor','k','LineWidth',W,'FaceAlpha',0.4);hold on;
text(1.95,y1-dy/2,{'under natural condition'},'Color','k','FontSize',14);
text(1.95,y2-dy/2,{'with mechanical compression'},'Color','k','FontSize',14);
DescriptorName={'General Sphericity','Diameter Sphericity','Intercept Sphericity','Maximum Projection Sphericity','Hayakawa Roundness','Spreading Index',...
    'Elongation Ratio','Pivotability Index','Wilson Flatness Index','Hayakawa Flatness Ratio','Huang Shape Factor','Corey Shape Factor'};
for DesNum=1:12
    plot([DesNum,DesNum+0.8],[0.265,0.265],'Color','k','LineWidth',W);hold on;
    text(DesNum+0.1,0.285,DescriptorName{1,DesNum},'Color','k','FontSize',14,'Rotation',35);
end

xticks([]);yticks(0:0.1:0.3);
ylabel({'Cell morphology variation ','at ~150 min after 4-cell stage'});
axis([0.5 12.8+0.5 0 0.4]);
set(gca,'FontSize',14,'FontName','Arial');
set(gca,'unit','centimeters','position',[3 1.5 14 8]);
set(gcf,'unit','centimeters','position',[2 8 20 12]);

