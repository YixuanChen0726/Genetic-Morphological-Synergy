% 画图，不同的命运分化对应的基因频率
clear all;clc;

% list the file name of all the gene expession data
filelist=dir(['..\..\bin\Gene Expression Data\Union Data\*.mat']);
genenamelist={};
for i=1:length(filelist)
    filename=filelist(i).name;
    namelist=strsplit(filename,{'_','.mat'});
    genenamelist=[genenamelist;namelist(1,1)];
end
genename=unique(genenamelist,'stable');

% set table
CellDiffTable=cell(22,length(genename)+1);
CellDiffTable(1,2:end)=genename';

% cell list
load(['..\..\bin\Coasymmetry\AsymCluster.mat']);
CellList=SortCellFreq(1:90,1);

% fate list
AllFateList={'Neuron';'Pharynx';'Muscle';'Skin';'Intestine';'Germ Cell'};

LineIndex=1;
for i=1:length(AllFateList)
    Fate1Name=AllFateList{i,1};
    for j=i:length(AllFateList)
        Fate2Name=AllFateList{j,1};
        TitleName=[Fate1Name,' / ',Fate2Name];
        LineIndex=LineIndex+1;
        CellDiffTable{LineIndex,1}=TitleName;
    end
end

% initial fequency
CellDiffTable(2:end,2:end)=num2cell(zeros(size(CellDiffTable,1)-1,size(CellDiffTable,2)-1));

for i=1:length(CellList)
    CellName=CellList{i,1};
    load(['..\..\bin\Fate Table\',CellName,'_FateTable.mat']);
    for j=2:size(CellDiffTable,1)
        string=split(CellDiffTable(j,1),' / ');
        Fate1=string{1,1};Fate2=string{2,1};
        % anterior 子代命运
        Row=find(strcmp(FateTable(:,1),Fate1));
        % posterior 子代命运
        Column=find(strcmp(FateTable(1,:),Fate2));
        if Row~=Column
            GeneList1=FateTable{Row,Column};
            GeneList2=FateTable{Column,Row};
            GeneList=[GeneList1;GeneList2];
        else
            GeneList=FateTable{Row,Column};
        end
        if ~isempty(GeneList)
            for k=1:length(GeneList)
                GeneName=GeneList{k,1};
                a=find(strcmp(CellDiffTable(1,:),GeneName));
                CellDiffTable{j,a}=CellDiffTable{j,a}+1;
            end
        end
    end
end
save(['Gene_Diff_Freq.mat'],'CellDiffTable','-v7.3');

