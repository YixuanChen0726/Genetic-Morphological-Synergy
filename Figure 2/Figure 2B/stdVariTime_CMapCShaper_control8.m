clear all;clc;close all;

% CMap
LifeCycleCMap=readcell(['..\..\bin\DataSet CMap\Life Cycle\Membrane\MeanEff_WT_Membrane_LifeCycle.csv']);
% CShaper
LifeCycleCShaper=readcell(['..\..\bin\DataSet CShaper\Membrane Life Cycle\MeanEff_LifeCycle.csv']);
LifeCycleCShaper(2:end,2:end)=num2cell(cell2mat(LifeCycleCShaper(2:end,2:end))*1.39);

% the last frame from 4-cell stage
ABaRow=find(strcmp(LifeCycleCMap(1,:),'ABa'));ABpRow=find(strcmp(LifeCycleCMap(1,:),'ABp'));
EMSRow=find(strcmp(LifeCycleCMap(1,:),'EMS'));P2Row=find(strcmp(LifeCycleCMap(1,:),'P2'));
StartTime=min([LifeCycleCMap{3,ABaRow},LifeCycleCMap{3,ABpRow},LifeCycleCMap{3,EMSRow},LifeCycleCMap{3,P2Row}]);
LifeCycleCMap(2:end,2:end)=num2cell(cell2mat(LifeCycleCMap(2:end,2:end))-StartTime);
LifeCycleCShaper(2:end,2:end)=num2cell(cell2mat(LifeCycleCShaper(2:end,2:end))-StartTime);
% shared cells in CMap and CShaper
CellList=intersect(LifeCycleCMap(1,2:end),LifeCycleCShaper(1,2:end));
LifeCycle=cell(3,1);LifeCycle(:,1)=LifeCycleCMap(:,1);
for CellIndex=1:length(CellList)
    Row=find(strcmp(LifeCycleCMap(1,:),CellList{CellIndex}));
    LifeCycle=[LifeCycle,LifeCycleCMap(:,Row)];
end

% General Sphericity
load(['..\..\bin\DataSet CMap\General Sphericity\DataRevisedScale\MeanRevise_Sphericity.mat']);
GSphStd=EffCell(meanGSph,LifeCycle);

% Diameter Sphericity
load(['v\Diameter Sphericity\DataRevisedScale\MeanRevise_DSphericity.mat']);
DSphStd=EffCell(meanDSph,LifeCycle);

% Intercept Sphericity
load(['..\..\bin\DataSet CMap\Intercept Sphericity\DataRevisedScale\MeanRevise_ISphericity.mat']);
ISphStd=EffCell(meanISph,LifeCycle);

% Maximum Projection Sphericity
load(['..\..\bin\DataSet CMap\Maximum Projection Sphericity\DataRevisedScale\MeanRevise_MPSphericity.mat']);
MPSphStd=EffCell(meanMPSph,LifeCycle);

% Hayakawa Roundness
load(['..\..\bin\DataSet CMap\Hayakawa Roundness\DataRevisedScale\MeanRevise_Roundness.mat']);
RoundStd=EffCell(meanRound,LifeCycle);

% Spreading Index
load(['..\..\bin\DataSet CMap\Spreading Index\DataRevisedScale\MeanRevise_SpreadingIndex.mat']);
SpreadStd=EffCell(meanSpread,LifeCycle);

% Elongation Ratio
load(['..\..\bin\DataSet CMap\Elongation Ratio\DataRevisedScale\MeanRevise_Elongation.mat']);
ElongationStd=EffCell(meanElongation,LifeCycle);

% Pivotability Index
load(['..\..\bin\DataSet CMap\Pivotability Index\DataRevisedScale\MeanRevise_Pivotability.mat']);
PivotabilityStd=EffCell(meanPivo,LifeCycle);

% Wilson Flatness Index
load(['..\..\bin\DataSet CMap\Wilson Flatness Index\DataRevisedScale\MeanRevise_WFlatness.mat']);
WFlatStd=EffCell(meanWFlat,LifeCycle);

% Hayakawa Flatness Ratio
load(['..\..\bin\DataSet CMap\Hayakawa Flatness Ratio\DataRevisedScale\MeanRevise_HFlatness.mat']);
HFlatStd=EffCell(meanHFlat,LifeCycle);

% Huang Shape Factor
load(['..\..\bin\DataSet CMap\Huang Shape Factor\DataRevisedScale\MeanRevise_Huang.mat']);
HuangStd=EffCell(meanHuang,LifeCycle);

% Corey Shape Factor
load(['..\..\bin\DataSet CMap\Corey Shape Factor\DataRevisedScale\MeanRevise_Corey.mat']);
CoreyStd=EffCell(meanCorey,LifeCycle);

% merge data
CellVari=cell(size(GSphStd,1)+1,14);
CellVari(:,1)=[{[]};GSphStd(:,1)];
CellVari(1,2:end)={'Time Point','General Sphericity','Diameter Sphericity','Intercept Sphericity','Maximum Projection Sphericity',...
    'Hayakawa Roundness','Spreading Index','Elongation Ratio','Pivotability Index',...
    'Wilson Flatness Index','Hayakawa Flatness Ratio','Huang Shape Factor','Corey Shape Factor'};
CellVari(2:end,2)=GSphStd(:,2);
CellVari(2:end,3:end)=[GSphStd(:,3),DSphStd(:,3),ISphStd(:,3),MPSphStd(:,3),RoundStd(:,3),SpreadStd(:,3),ElongationStd(:,3),PivotabilityStd(:,3),WFlatStd(:,3),HFlatStd(:,3),HuangStd(:,3),CoreyStd(:,3)];

save(['CellVariationTime_CMapCShaper_control8.mat'],'CellVari','-v7.3');

% screen effective cell
function Modifiedstd=EffCell(Descriptor,LifeCycle)
% cell name, time ,std
Modifiedstd=cell(1,3);LinIndex=0;
for CellNum=2:size(LifeCycle,2)
    CellName=LifeCycle{1,CellNum};
    CellNameRow=find(strcmp(Descriptor(:,1),CellName));
    CellStd=cell2mat(Descriptor(CellNameRow,5));
    CellBirth=LifeCycle{2,CellNum};CellDeath=LifeCycle{3,CellNum};
    TimeList=linspace(CellBirth,CellDeath,length(CellNameRow));
    for Time=1:length(CellNameRow)
        LinIndex=LinIndex+1;
        Modifiedstd(LinIndex,:)={CellName,TimeList(Time),CellStd(Time)};
    end
end
end