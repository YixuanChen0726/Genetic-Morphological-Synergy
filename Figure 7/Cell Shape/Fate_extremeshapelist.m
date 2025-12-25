% For each fate, select the top five from each indicator.
clear all;clc;close all;

% the fate table
load(['..\..\bin\FateClassification.mat']);

FateNameList={'Neuron','Pharynx','Muscle','Skin','Intestine'};

for i=1:length(FateNameList)
    % Extreme Fate List
    ExtremeFateList=cell(3,12);
    ExtremeFateList={'general sphericity','diameter sphericity','intercept sphericity','maximum projection sphericity',...
        'Hayakawa roundness','spreading index','elongation ratio','pivotability index',...
        'Wilson flatness index','Hayakawa flatness ratio','Huang shape factor','Corey shape factor'};
    FateName=FateNameList{1,i};
    % general sphericity
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\General Sphericity\DataCleanLife','Sphericity');
    ExtremeFateList([2,3],1)={HighData;LowData};
    % diameter sphericity
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Diameter Sphericity\DataCleanLife','DSphericity');
    ExtremeFateList([2,3],2)={HighData;LowData};
    % intercept sphericity
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Intercept Sphericity\DataCleanLife','ISphericity');
    ExtremeFateList([2,3],3)={HighData;LowData};
    % maximum projection sphericity
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Maximum Projection Sphericity\DataCleanLife','MPSphericity');
    ExtremeFateList([2,3],4)={HighData;LowData};
    % Hayakawa roundness
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Hayakawa Roundness\DataCleanLife','HayakawaRoundness');
    ExtremeFateList([2,3],5)={HighData;LowData};
    % spreading index
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Spreading Index\DataCleanLife','SpreadingIndex');
    ExtremeFateList([2,3],6)={HighData;LowData};
    % elongation ratio
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Elongation Ratio\DataCleanLife','Elongation');
    ExtremeFateList([2,3],7)={HighData;LowData};
    % pivotability index
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Pivotability Index\DataCleanLife','Pivotability');
    ExtremeFateList([2,3],8)={HighData;LowData};
    % Wilson flatness index
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Wilson Flatness Index\DataCleanLife','Flatness');
    ExtremeFateList([2,3],9)={HighData;LowData};
    % Hayakawa flatness ratio
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Hayakawa Flatness Ratio\DataCleanLife','Flatness');
    ExtremeFateList([2,3],10)={HighData;LowData};
    % Huang shape factor
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Huang Shape Factor\DataCleanLife','Huang');
    ExtremeFateList([2,3],11)={HighData;LowData};
    % Corey shape factor
    [HighData,LowData]=ExtremeList(FateClas,FateName,'..\..\bin\DataSet CMap\Corey Shape Factor\DataCleanLife','Corey');
    ExtremeFateList([2,3],12)={HighData;LowData};
    save([FateName,'_ExtremeShapeList.mat'],'ExtremeFateList','-v7.3');
    if strcmp(FateName,'Skin')
        save(['Ectoderm_ExtremeShapeList.mat'],'ExtremeFateList','-v7.3');
    end
    if strcmp(FateName,'Intestine')
        save(['Entoderm_ExtremeShapeList.mat'],'ExtremeFateList','-v7.3');
    end
end

% shape descriptor
function [HighData,LowData]=ExtremeList(FateClas,FateName,FileName,ShapeName)
AllShapeData=[];
for SampleNum=1:8
    val_struct=load([FileName,'\WT_Sample',num2str(SampleNum),'_',ShapeName,'.mat']);
    val_names=fieldnames(val_struct);
    ShapeData=getfield(val_struct,val_names{1});
    % [SampleNum,Cell,Frame,ShapeData]
    AllShapeData=[AllShapeData;ShapeData(:,[1,2,4,5])];
end

% fate list
CellList=FateClas(find(strcmp(FateClas(:,2),FateName)),1);

FateShapeList=[];
for i=1:length(CellList)
    Row=find(strcmp(AllShapeData(:,2),CellList{i,1}));
    FateShapeList=[FateShapeList;AllShapeData(Row,:)];
end

[B,I]=sort(cell2mat(FateShapeList(:,4)));
SortShapeData=FateShapeList(I,:);

LowData=SortShapeData(1:5,:);
HighData=SortShapeData(end-4:end,:);
end
