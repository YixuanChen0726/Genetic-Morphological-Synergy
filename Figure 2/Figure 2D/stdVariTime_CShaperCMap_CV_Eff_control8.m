% morphology descriptor/ cell list
clear all;clc;

% WT,CMap
LifeCycleCMap=readcell(['..\..\bin\DataSet CMap\Life Cycle\Membrane\MeanEff_WT_Membrane_LifeCycle.csv']);
% CShaper
LifeCycleCShaper=readcell(['..\..\bin\DataSet CShaper\Membrane Life Cycle\MeanEff_LifeCycle.csv']);
LifeCycleCShaper(2:end,2:end)=num2cell(cell2mat(LifeCycleCShaper(2:end,2:end))*1.39);

% the last frame of 4-cell stage
ABaRow=find(strcmp(LifeCycleCShaper(1,:),'ABa'));ABpRow=find(strcmp(LifeCycleCShaper(1,:),'ABp'));
EMSRow=find(strcmp(LifeCycleCShaper(1,:),'EMS'));P2Row=find(strcmp(LifeCycleCShaper(1,:),'P2'));
StartTime=min([LifeCycleCShaper{3,ABaRow},LifeCycleCShaper{3,ABpRow},LifeCycleCShaper{3,EMSRow},LifeCycleCShaper{3,P2Row}]);
LifeCycleCShaper(2:end,2:end)=num2cell(cell2mat(LifeCycleCShaper(2:end,2:end))-StartTime);
% shared cells of CMap and CShaper
CellList=intersect(LifeCycleCMap(1,2:end),LifeCycleCShaper(1,2:end));
LifeCycle=cell(3,1);LifeCycle(:,1)=LifeCycleCMap(:,1);
for CellIndex=1:length(CellList)
    Row=find(strcmp(LifeCycleCShaper(1,:),CellList{CellIndex}));
    LifeCycle=[LifeCycle,LifeCycleCShaper(:,Row)];
end

% CV
load(['..\..\bin\DataSet CShaper\Cell Contact\DataRevisedScale\Mean_CV_CShaper_Eff_control8.mat']);
CVStd=EffCell(CV,LifeCycle);

% merge data
CellVari=cell(size(CVStd,1)+1,3);
CellVari(:,1)=[{[]};CVStd(:,1)];
CellVari(1,2:end)={'Time Point','CV'};
CellVari(2:end,2)=CVStd(:,2);
CellVari(2:end,3)=[CVStd(:,3)];

save(['CellVariationTime_CShaperCMap_ContactCV_control8.mat'],'CellVari','-v7.3');

% screen effective cells
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