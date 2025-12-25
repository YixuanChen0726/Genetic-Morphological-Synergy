% morphology descriptor/ time
clear all;clc;

% CShaper
LifeCycleCShaper=readcell(['..\..\bin\DataSet CShaper\Membrane Life Cycle\MeanEff_LifeCycle.csv']);
LifeCycleCShaper(2:end,2:end)=num2cell(cell2mat(LifeCycleCShaper(2:end,2:end))*1.39);
% CMap
LifeCycleCMap=readcell(['..\..\bin\DataSet CMap\Life Cycle\Membrane\MeanEff_WT_Membrane_LifeCycle.csv']);

% the last frame of 4-cell stage
ABaRow=find(strcmp(LifeCycleCShaper(1,:),'ABa'));ABpRow=find(strcmp(LifeCycleCShaper(1,:),'ABp'));
EMSRow=find(strcmp(LifeCycleCShaper(1,:),'EMS'));P2Row=find(strcmp(LifeCycleCShaper(1,:),'P2'));
StartTime=min([LifeCycleCShaper{3,ABaRow},LifeCycleCShaper{3,ABpRow},LifeCycleCShaper{3,EMSRow},LifeCycleCShaper{3,P2Row}]);
LifeCycleCShaper(2:end,2:end)=num2cell(cell2mat(LifeCycleCShaper(2:end,2:end))-StartTime);
LifeCycleCMap(2:end,2:end)=num2cell(cell2mat(LifeCycleCMap(2:end,2:end))-StartTime);
% cells shared in CMap and CShaper
CellList=intersect(LifeCycleCShaper(1,2:end),LifeCycleCMap(1,2:end));
LifeCycle=cell(3,1);LifeCycle(:,1)=LifeCycleCShaper(:,1);
for CellIndex=1:length(CellList)
    Row=find(strcmp(LifeCycleCShaper(1,:),CellList{CellIndex}));
    LifeCycle=[LifeCycle,LifeCycleCShaper(:,Row)];
end

% General Sphericity
load(['..\..\bin\DataSet CShaper\General Sphericity\MeanRevise_Sphericity_control8.mat']);
GSphStd=EffCell(meanGSph,LifeCycle);

% Diameter Sphericity
load(['..\..\bin\DataSet CShaper\Diameter Sphericity\MeanRevise_DSphericity_control8.mat']);
DSphStd=EffCell(meanDSph,LifeCycle);

% Intercept Sphericity
load(['..\..\bin\DataSet CShaper\Intercept Sphericity\MeanRevise_ISphericity_control8.mat']);
ISphStd=EffCell(meanISph,LifeCycle);

% Maximum Projection Sphericity
load(['..\..\bin\DataSet CShaper\Maximum Projection Sphericity\MeanRevise_MPSphericity_control8.mat']);
MPSphStd=EffCell(meanMPSph,LifeCycle);

% Hayakawa Roundness
load(['..\..\bin\DataSet CShaper\Hayakawa Roundness\MeanRevise_Roundness_control8.mat']);
RoundStd=EffCell(meanRound,LifeCycle);

% Spreading Index
load(['..\..\bin\DataSet CShaper\Spreading Index\MeanRevise_SpreadingIndex_control8.mat']);
SpreadStd=EffCell(meanSpread,LifeCycle);

% Elongation Ratio
load(['..\..\bin\DataSet CShaper\Elongation Ratio\MeanRevise_Elongation_control8.mat']);
ElongationStd=EffCell(meanElongation,LifeCycle);

% Pivotability Index
load(['..\..\bin\DataSet CShaper\Pivotability Index\MeanRevise_Pivotability_control8.mat']);
PivotabilityStd=EffCell(meanPivotability,LifeCycle);

% Wilson Flatness Index
load(['..\..\bin\DataSet CShaper\Wilson Flatness Index\MeanRevise_WFlatness_control8.mat']);
WFlatStd=EffCell(meanWFlat,LifeCycle);

% Hayakawa Flatness Ratio
load(['..\..\bin\DataSet CShaper\Hayakawa Flatness Ratio\MeanRevise_HFlatness_control8.mat']);
HFlatStd=EffCell(meanHFlat,LifeCycle);

% Huang Shape Factor
load(['..\..\bin\DataSet CShaper\Huang Shape Factor\MeanRevise_Huang_control8.mat']);
HuangStd=EffCell(meanHuang,LifeCycle);

% Corey Shape Factor
load(['..\..\bin\DataSet CShaper\Corey Shape Factor\MeanRevise_Corey_control8.mat']);
CoreyStd=EffCell(meanCorey,LifeCycle);

% merge data
CellVari=cell(size(GSphStd,1)+1,14);
CellVari(:,1)=[{[]};GSphStd(:,1)];
CellVari(1,2:end)={'Time Point','General Sphericity','Diameter Sphericity','Intercept Sphericity','Maximum Projection Sphericity',...
    'Hayakawa Roundness','Spreading Index','Elongation Ratio','Pivotability Index',...
    'Wilson Flatness Index','Hayakawa Flatness Ratio','Huang Shape Factor','Corey Shape Factor'};
CellVari(2:end,2)=GSphStd(:,2);
CellVari(2:end,3:end)=[GSphStd(:,3),DSphStd(:,3),ISphStd(:,3),MPSphStd(:,3),RoundStd(:,3),SpreadStd(:,3),ElongationStd(:,3),PivotabilityStd(:,3),WFlatStd(:,3),HFlatStd(:,3),HuangStd(:,3),CoreyStd(:,3)];

save(['CellVariationTime_CShaperCMap_control8.mat'],'CellVari','-v7.3');

% screen data
function Modifiedstd=EffCell(Descriptor,LifeCycle)
% cell name, time, std
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