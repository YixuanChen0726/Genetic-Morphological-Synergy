% Statistical differentiation event
clear all; clc; close all;

% five fates:Neuron,Pharynx,Muscle,Skin,Intestine
% For each candidate cell, determine the composition of its fate, and compare the differentiation of all progeny cells with the differentiation of individual fates.
% If the differentiation differences in individual fates become greater, the average degree of asymmetry also increases.
% This gene can be regarded as a marker for a certain fate.

% cell list
load(['..\..\bin\Coasymmetry\AsymCluster.mat']);
CellList=SortCellFreq(1:90,1);

% lineage tree
load(['..\..\bin\LineageTree.mat']);
% fate list
AllFateList={'Neuron';'Pharynx';'Muscle';'Skin';'Intestine'};
% fate distribution
load('..\..\bin\FateClassification.mat');

for i=1:length(CellList)
    % table of different fates
    FateTable=cell(2,length(AllFateList));
    FateTable(1,:)=AllFateList';
    FateTable(2,:)={{},{},{},{},{}};
    % mother cell
    CellName=CellList{i,1};
    % two sublineage daughters
    [a,b]=find(strcmp(LineageTree,CellName));
    CellLayer=b+1;Cell1Index=a*2-1;Cell2Index=a*2;
    Dau1Name=LineageTree{Cell1Index,CellLayer};
    Dau2Name=LineageTree{Cell2Index,CellLayer};
    % Daughter1/2 cell list
    Dau1List={};Dau2List={};
    for DauLayer=CellLayer:1:size(LineageTree,2)
        n=2^(DauLayer-CellLayer);
        for Dau1Index=(n*Cell1Index-(n-1)):1:min([(n*Cell1Index),size(LineageTree,1)])
            DauName=LineageTree{Dau1Index,DauLayer};
            if all(DauName)
                Dau1List=[Dau1List;{DauName}];
            end
        end
        for Dau2Index=(n*Cell2Index-(n-1)):1:min([(n*Cell2Index),size(LineageTree,1)])
            DauName=LineageTree{Dau2Index,DauLayer};
            if all(DauName)
                Dau2List=[Dau2List;{DauName}];
            end
        end
    end
    % all the fates of two daughters
    [data,a,b]=intersect(FateClas(:,1),Dau1List);
    Dau1AllFate=unique(FateClas(a,2));
    Dau1Fate=intersect(Dau1AllFate,AllFateList);
    [data,a,b]=intersect(FateClas(:,1),Dau2List);
    Dau2AllFate=unique(FateClas(a,2));
    Dau2Fate=intersect(Dau2AllFate,AllFateList);

    % [gene name, p value, asymmetry]
    load(['..\..\bin\Coasymmetry\all gene\',CellName,'_AllGeneSig.mat']);
    AllGeneData=AllGeneSigni([1,5,6],:);

    % the difference of specific fates
    for Fate1=1:length(Dau1Fate)
        DFate1=Dau1Fate{Fate1,1};
        for Fate2=1:length(Dau2Fate)
            DFate2=Dau2Fate{Fate2,1};
            if ~strcmp(DFate1,DFate2)
                load(['..\..\bin\Coasymmetry\partial fate gene\',CellName,'_',DFate1,'_',DFate2,'.mat']);
                ParGeneData=PartialGeneSigni([1,5,6],:);
                % 如果满足增强条件，放置在FateTable相应位置
                for GeneNum=1:size(AllGeneData,2)
                    GeneName=AllGeneData{1,GeneNum};
                    AllPValue=AllGeneData{2,GeneNum};
                    AllAsym=AllGeneData{3,GeneNum};
                    ParPValue=ParGeneData{2,GeneNum};
                    ParAsym=ParGeneData{3,GeneNum};
                    % 双强
                    if (ParPValue<AllPValue) && (abs(AllAsym)<abs(ParAsym))
                        if ParAsym>0
                            n=find(strcmp(FateTable(1,:),DFate1));
                            FateTable{2,n}=[FateTable{2,n};GeneName];
                        elseif ParAsym<0
                            n=find(strcmp(FateTable(1,:),DFate2));
                            FateTable{2,n}=[FateTable{2,n};GeneName];
                        end
                    end
                end
            end
        end
    end
    save(['Statistical FateMarker\',CellName,'_FateGeneMarker.mat'],'FateTable','-v7.3');
end
