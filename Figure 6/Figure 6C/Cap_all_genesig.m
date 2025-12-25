% Cap sublineage differentiation
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

% list the file name of all data
filelist=dir(['..\..\bin\Gene Expression Data\Union Data\*.mat']);
genenamelist={};
for i=1:length(filelist)
    filename=filelist(i).name;
    namelist=strsplit(filename,{'_','.mat'});
    genenamelist=[genenamelist;namelist(1,1)];
end
genename=unique(genenamelist,'stable');

AnteriorList=CapaList;PosteriorList=CappList;

% [gene name, anterior cells' expression, posterior cells' expression, mean expression of each sublineage, p value, asymmetry]
GeneSigni=cell(6,length(genename));

for i=1:length(genename)
    GeneSigni{1,i}=genename{i,1};
    Row1=find(strcmp(genenamelist,genename{i,1}));
    genedata=ArrData(filelist,Row1);
    Temp1={};
    for j=1:length(AnteriorList)
        CellName=AnteriorList{j,1};
        c=find(strcmp(genedata(:,1),CellName));
        Temp1=[Temp1;genedata(c,[1,3])];
    end
    GeneSigni{2,i}=Temp1;
    Temp2={};
    for j=1:length(PosteriorList)
        CellName=PosteriorList{j,1};
        c=find(strcmp(genedata(:,1),CellName));
        Temp2=[Temp2;genedata(c,[1,3])];
    end
    GeneSigni{3,i}=Temp2;
    % delete NaN
    ExpData1=cell2mat(GeneSigni{2,i}(:,2));
    ExpData2=cell2mat(GeneSigni{3,i}(:,2));
    Row1=find(isnan(ExpData1));
    Dau1Data=ExpData1;Dau1Data(Row1,:)=[];
    Row2=find(isnan(ExpData2));
    Dau2Data=ExpData2;Dau2Data(Row2,:)=[];
    GeneSigni{4,i}=[mean(Dau1Data);mean(Dau2Data)]; 
    % p value and asymmetry
    [h,p,ci,stats]=ttest2(Dau1Data,Dau2Data);
    GeneSigni{5,i}=p;
    meandata1=GeneSigni{4,i}(1,1);
    meandata2=GeneSigni{4,i}(2,1);
    if any([meandata1,meandata2])
        asym=(meandata1-meandata2)/(meandata1+meandata2);
    else
        asym=nan;
    end
    GeneSigni{6,i}=asym;
    disp(i);
end
save(['Cap_all_genesig.mat'],'GeneSigni','-v7.3');

% p_threshold < 0.1
function genedata=ArrData(filelist,Row)
% data storage:[cell name, time, expression data]
filename=filelist(Row(1)).name;
load(['E:\Project-Gene Expression\Gene Expression Data\plotexpression\AllData\',filename]);
AllData=MeanExp(:,[1,2,4]);
if length(Row)>1
    for r=2:length(Row)
        filename=filelist(Row(r)).name;
        % gene table
        load(['E:\Project-Gene Expression\Gene Expression Data\plotexpression\AllData\',filename]);
        AllData=[AllData,MeanExp(:,4)];
    end
end
allexp=cell2mat(AllData(:,3:end));
Exp=mean(allexp,2,'omitnan');Exp(Exp<0.1)=0;
genedata=[AllData(:,[1,2]),num2cell(Exp)];
end
