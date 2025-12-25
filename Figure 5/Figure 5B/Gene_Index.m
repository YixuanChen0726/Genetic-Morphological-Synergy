% Gene index
clear all;clc; close all;

% Num2: Unique notation for a specific labeled gene.
NameDic=readcell(['..\..\bin\Gene Expression Data\Name Dictionary.csv']);
GeneList=unique(NameDic(2:end,2));

GeneDic=cell(0,2);LineIndex=0;
for GeneNum=1:length(GeneList)
    GeneName=GeneList{GeneNum,1};
    Row=find(strcmp(NameDic(:,2),GeneName));
    FileName=NameDic{Row(1),1};
    NumList=strsplit(FileName,{'_','.'});
    GeneIndex=str2double(NumList{1,3});
    LineIndex=LineIndex+1;
    GeneDic(LineIndex,:)={GeneName,GeneIndex};
end

writecell(GeneDic,['GeneIndex_Dictionary.csv']);