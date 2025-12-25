% For each fate, select the top ten from each indicator.
clear all;clc;close all;

% fate table
load(['E:\Project-Gene Expression\Organ Shape\FateClassification.mat']);

% Extreme Fate List
ExtremeFateList=cell(3,12);
ExtremeFateList={'general sphericity','diameter sphericity','intercept sphericity','maximum projection sphericity',...
    'Hayakawa roundness','spreading index','elongation ratio','pivotability index',...
    'Wilson flatness index','Hayakawa flatness ratio','Huang shape factor','Corey shape factor'};

% general sphericity
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\General Sphericity\DataCleanLife','Sphericity');
ExtremeFateList([2,3],1)={HighData;LowData};
% diameter sphericity
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Diameter Sphericity\DataCleanLife','DSphericity');
ExtremeFateList([2,3],2)={HighData;LowData};
% intercept sphericity
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Intercept Sphericity\DataCleanLife','ISphericity');
ExtremeFateList([2,3],3)={HighData;LowData};
% maximum projection sphericity
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Maximum Projection Sphericity\DataCleanLife','MPSphericity');
ExtremeFateList([2,3],4)={HighData;LowData};
% Hayakawa roundness
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Hayakawa Roundness\DataCleanLife','HayakawaRoundness');
ExtremeFateList([2,3],5)={HighData;LowData};
% spreading index
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Spreading Index\DataCleanLife','SpreadingIndex');
ExtremeFateList([2,3],6)={HighData;LowData};
% elongation ratio
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Elongation Ratio\DataCleanLife','Elongation');
ExtremeFateList([2,3],7)={HighData;LowData};
% pivotability index
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Pivotability Index\DataCleanLife','Pivotability');
ExtremeFateList([2,3],8)={HighData;LowData};
% Wilson flatness index
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Wilson Flatness Index\DataCleanLife','Flatness');
ExtremeFateList([2,3],9)={HighData;LowData};
% Hayakawa flatness ratio
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Hayakawa Flatness Ratio\DataCleanLife','Flatness');
ExtremeFateList([2,3],10)={HighData;LowData};
% Huang shape factor
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Huang Shape Factor\DataCleanLife','Huang');
ExtremeFateList([2,3],11)={HighData;LowData};
% Corey shape factor
[HighData,LowData]=ExtremeList(FateClas,'..\..\bin\DataSet CMap\Corey Shape Factor\DataCleanLife','Corey');
ExtremeFateList([2,3],12)={HighData;LowData};
save(['Mesoderm_ExtremeShapeList.mat'],'ExtremeFateList','-v7.3');

% shape descriptor
function [HighData,LowData]=ExtremeList(FateClas,FileName,ShapeName)
AllShapeData=[];
for SampleNum=1:8
    val_struct=load([FileName,'\WT_Sample',num2str(SampleNum),'_',ShapeName,'.mat']);
    val_names=fieldnames(val_struct);
    ShapeData=getfield(val_struct,val_names{1});
    % 1,2,4,5列，SampleNum，Cell，Frame，ShapeData
    AllShapeData=[AllShapeData;ShapeData(:,[1,2,4,5])];
end

CellList1=FateClas(find(strcmp(FateClas(:,2),'Neuron')),1);
CellList2=FateClas(find(strcmp(FateClas(:,2),'Pharynx')),1);
CellList3=FateClas(find(strcmp(FateClas(:,2),'Muscle')),1);
CellList=[CellList1;CellList2;CellList3];

FateShapeList=[];
for i=1:length(CellList)
    Row=find(strcmp(AllShapeData(:,2),CellList{i,1}));
    FateShapeList=[FateShapeList;AllShapeData(Row,:)];
end

[B,I]=sort(cell2mat(FateShapeList(:,4)));
SortShapeData=FateShapeList(I,:);

LowData=SortShapeData(1:10,:);
HighData=SortShapeData(end-9:end,:);
end
