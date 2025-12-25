% Integrate all the genetic markers related to fate
clear all;clc;close all;

GeneMarkerTab=cell(3,5);
GeneMarkerTab={'Neuron','Pharynx','Muscle','Skin','Intestine'};

filelist=dir(['Statistical FateMarker\*.mat']);
NeuronMarker=[];PharynxMarker=[];MuscleMarker=[];SkinMarker=[];IntestineMarker=[];
for i=1:length(filelist)
    filename=filelist(i).name;
    load(['Statistical FateMarker\',filename]);
    NeuronMarker=[NeuronMarker;FateTable{2,1}];
    PharynxMarker=[PharynxMarker;FateTable{2,2}];
    MuscleMarker=[MuscleMarker;FateTable{2,3}];
    SkinMarker=[SkinMarker;FateTable{2,4}];
    IntestineMarker=[IntestineMarker;FateTable{2,5}];
end

GeneMarkerTab{2,1}=unique(NeuronMarker);
GeneMarkerTab{2,2}=unique(PharynxMarker);
GeneMarkerTab{2,3}=unique(MuscleMarker);
GeneMarkerTab{2,4}=unique(SkinMarker);
GeneMarkerTab{2,5}=unique(IntestineMarker);

% Gene Dictionary
GeneDic=readcell('GeneIndex_Dictionary.csv');
% transfer to index
for i=1:size(GeneMarkerTab,2)
    GeneList=GeneMarkerTab{2,i};
    GeneIndexList=[];
    for j=1:length(GeneList)
        GeneName=GeneList{j,1};
        k=find(strcmp(GeneDic(:,1),GeneName));
        GeneIndex=GeneDic{k,2};
        GeneIndexList=[GeneIndexList;GeneIndex];
    end
    GeneMarkerTab{3,i}=GeneIndexList;
end
save(['Fate_GeneMarker.mat'],'GeneMarkerTab','-v7.3');
