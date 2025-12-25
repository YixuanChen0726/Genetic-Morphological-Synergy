% Conduct significance analysis on all shapes of all cells
clear all;clc;

% lineage tree
load(['..\..\bin\LineageTree.mat']);
AllCellName={};
for i=1:size(LineageTree,1)
    for j=1:size(LineageTree,2)
        CellName=LineageTree{i,j};
        if all(CellName)
            AllCellName=[AllCellName;{CellName}];
        end
    end
end
% Capa: all cells
CapaList=AllCellName(find(contains(AllCellName,'Capa')),1);
% Capp: all cells
CappList=AllCellName(find(contains(AllCellName,'Capp')),1);

% raw data
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\GSphTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\DSphTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\ISphTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\MPSphTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\HRoundTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\SpreadTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\ElongationTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\SpreadTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\PivotabilityTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\HFlatTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\WFlatTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\HuangTree_Shape.mat']);
load(['..\..\bin\Coasymmetry\Sublineage Asymmetry\sister sublineage shape\CoreyTree_Shape.mat']);

% [gene name, anterior cells' expression, posterior cells' expression, mean expression of each sublineage, p value, asymmetry]
ShapeSigni=cell(6,12);
ShapeSigni(1,:)={'general sphericity','diameter sphericity','intercept sphericity','maximum projection sphericity',...
    'Hayakawa roundness','spreading index','elongation ratio','pivotability index',...
    'Hayakawa Flatness ratio','Wilson Flatness index','Huang shape factor','Corey shape factor'};

[Temp1,Temp2,meanData,p,asym]=SigniFun(GSphTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,1)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(DSphTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,2)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(ISphTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,3)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(MPSphTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,4)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(HRoundTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,5)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(SpreadTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,6)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(ElonTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,7)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(PivoTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,8)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(HFlatTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,9)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(WFlatTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,10)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(HuangTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,11)={Temp1;Temp2;meanData;p;asym};
[Temp1,Temp2,meanData,p,asym]=SigniFun(CoreyTree,LineageTree,CapaList,CappList);
ShapeSigni(2:end,12)={Temp1;Temp2;meanData;p;asym};

save(['Cap_all_shapesig.mat'],'ShapeSigni','-v7.3');

function [Temp1,Temp2,meanData,p,asym]=SigniFun(DescriptorTree,LineageTree,AnteriorList,PosteriorList)
[a1,b1]=find(strcmp(LineageTree,'Capa'));
[a2,b2]=find(strcmp(LineageTree,'Capp'));
anData=DescriptorTree{a1,b1};poData=DescriptorTree{a2,b2};
Temp1={};
for j=1:length(AnteriorList)
    CellName=AnteriorList{j,1};
    c=find(strcmp(anData(:,1),CellName));
    Temp1=[Temp1;anData(c,:)];
end
Temp2={};
for j=1:length(PosteriorList)
    CellName=PosteriorList{j,1};
    c=find(strcmp(poData(:,1),CellName));
    Temp2=[Temp2;poData(c,:)];
end
DauData1=cell2mat(Temp1(:,2));DauData2=cell2mat(Temp2(:,2));
% mean
meandata1=mean(DauData1);meandata2=mean(DauData2);
meanData=[meandata1;meandata2];
% p value & asymmetry
[h,p,ci,stats]=ttest2(DauData1,DauData2);
asym=(meandata1-meandata2)/(meandata1+meandata2);
end
