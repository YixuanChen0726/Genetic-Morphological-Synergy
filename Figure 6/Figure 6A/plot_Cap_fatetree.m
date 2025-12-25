% plot lineage tree colored by fate
clear all;clc;close all;

% lineage tree
load(['..\..\bin\LineageTree.mat']);
% fate tree
load(['..\..\bin\FateClassification.mat']);
% Position Tree
load(['..\..\bin\CellPositionTree_All.mat']);
% cell life cycle
AllLifeCycle=readcell(['..\..\bin\DataSet CMap\Life Cycle\Nucleus\MeanTot_WT_Nucleus_LifeCycle.csv']);
% the last frame of 4-cell stage
ABaRow=find(strcmp(AllLifeCycle(1,:),'ABa'));ABpRow=find(strcmp(AllLifeCycle(1,:),'ABp'));
EMSRow=find(strcmp(AllLifeCycle(1,:),'EMS'));P2Row=find(strcmp(AllLifeCycle(1,:),'P2'));
StartTime=min([AllLifeCycle{3,ABaRow},AllLifeCycle{3,ABpRow},AllLifeCycle{3,EMSRow},AllLifeCycle{3,P2Row}]);
AllLifeCycle(2:end,2:end)=num2cell(cell2mat(AllLifeCycle(2:end,2:end))-StartTime);

% fate color
NeuronColor=[116,190,200]./255;
PharynxColor=[132,62,160]./255;
MuscleColor=[231,121,184]./255;
IntestineColor=[235,194,61]./255;
SkinColor=[161,210,88]./255;
GermColor=[136,169,228]./255;
DeathColor=[207,41,41]./255;
UnspecifiedColor=[185,185,185]./255;

% Cap and its offspring
k=find(contains(FateClas(:,1),'Cap'));
CellList=FateClas(k,1);

LineageData=cell(size(LineageTree));
for i=1:length(CellList)
    CellName=CellList{i,1};
    % cell fate
    Row=find(strcmp(FateClas(:,1),CellName));
    CellFate=FateClas{Row,2};
    % x position
    [a,b]=find(strcmp(LineageTree,CellName));
    CellPosition=PositionTree{a,b};
    % cell color
    switch CellFate
        case 'Neuron'
            CellColor=NeuronColor;
        case 'Pharynx'
            CellColor=PharynxColor;
        case 'Muscle'
            CellColor=MuscleColor;
        case 'Intestine'
            CellColor=IntestineColor;
        case 'Skin'
            CellColor=SkinColor;
        otherwise
            CellColor=UnspecifiedColor;
    end
    % cell life
    MotherRow=find(strcmp(AllLifeCycle(1,:),CellName));
    BirthTime=AllLifeCycle{2,MotherRow};
    Dau1Name=LineageTree{a*2-1,b+1};
    Dau2Name=LineageTree{a*2,b+1};
    Dau1Row=find(strcmp(AllLifeCycle(1,:),Dau1Name));
    Dau2Row=find(strcmp(AllLifeCycle(1,:),Dau2Name));
    if isempty(Dau1Row) && isempty(Dau2Row)
        DeathTime=AllLifeCycle{3,MotherRow};
    elseif isempty(Dau1Row)
        DeathTime=AllLifeCycle{2,Dau2Row};
    elseif isempty(Dau2Row)
        DeathTime=AllLifeCycle{2,Dau1Row};
    else
        Dau1Birth=AllLifeCycle{2,Dau1Row};
        Dau2Birth=AllLifeCycle{2,Dau2Row};
        DeathTime=max([Dau1Birth,Dau2Birth]);
    end
    % [x, time point, mean data, std data]
    CharaLineage={CellPosition,[BirthTime;DeathTime],CellColor};
    LineageData{a,b}=CharaLineage;
end

% plot
LinW=1.5;
for CellLayer=size(LineageData,2):-1:1
    for CellIndex=1:1:size(LineageData,1)
        Data=LineageData{CellIndex,CellLayer};
        if ~isempty(Data)
            Xp=Data{1,1};
            BirthTime=Data{1,2}(1,1);
            DeathTime=Data{1,2}(2,1);
            Color=Data{1,3};
            plot([Xp,Xp],[BirthTime,DeathTime],'Color',Color,'LineWidth',LinW);hold on;
            MData=LineageData{ceil(CellIndex/2),CellLayer-1};
            if ~isempty(MData)
                MXp=MData{1,1};
                plot([Xp,MXp],[BirthTime,BirthTime],'Color',Color,'LineWidth',LinW);hold on;
            end
        end
    end
end

% axis range
xrange=[];
for CellIndex=1:size(LineageData,1)
    for CellLayer=1:size(LineageData,2)
        Data=LineageData{CellIndex,CellLayer};
        if ~isempty(Data)
            Position=Data{1,1};
            xrange=[xrange;Position];
        end
    end
end

% mark Lineage
[Row,Column]=find(strcmp(LineageTree,'Cap'));
X=LineageData{Row,Column}{1,1};
Y=LineageData{Row,Column}{1,2}(1,1);
text(X,Y-18,'Cap','FontSize',15,'FontName','Arial','HorizontalAlignment','center');

% cell position
dau1name='Capa';dau2name='Capp';
[D1Row,D1Column]=find(strcmp(LineageTree,dau1name));
Dau1Position=[];
for DauLayer=D1Column:1:size(LineageTree,2)
    n=2^(DauLayer-D1Column);
    for DauIndex=(n*D1Row-(n-1)):1:min([(n*D1Row),size(LineageTree,1)])
        DauPos=PositionTree{DauIndex,DauLayer};
        Dau1Position=[Dau1Position;DauPos];
    end
end
Dau1range=[min(Dau1Position);max(Dau1Position)];
[D2Row,D2Column]=find(strcmp(LineageTree,dau2name));
Dau2Position=[];
for DauLayer=D2Column:1:size(LineageTree,2)
    n=2^(DauLayer-D2Column);
    for DauIndex=(n*D2Row-(n-1)):1:min([(n*D2Row),size(LineageTree,1)])
        DauPos=PositionTree{DauIndex,DauLayer};
        Dau2Position=[Dau2Position;DauPos];
    end
end
Dau2range=[min(Dau2Position);max(Dau2Position)];
% plot line
plot([Dau1range(1,1),Dau1range(2,1)],[300,300],'-','LineWidth',LinW,'Color','k');hold on;
plot([Dau2range(1,1),Dau2range(2,1)],[300,300],'-','LineWidth',LinW,'Color','k');hold on;

text(mean(Dau1range),310,{[num2str(6/8*100,'%.1f'),'%'],'ventral'},'FontSize',15,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','top');
text(mean(Dau2range),310,{[num2str(7/8*100,'%.1f'),'%'],'dorsal'},'FontSize',15,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','top');

axis([min(xrange)-2,max(xrange)+2,30,320]);
yticks(0:40:280);box off;
ylabel({'\it t\rm (min)'},'Rotation',0,'Position',[min(xrange)-6.5 330]);
set(gca,'YDir','reverse','YAxisLocation','left','ycolor','k');
set(gca,'xtick',[],'XAxisLocation','top','xcolor','none');
set(gca,'LineWidth',1,'FontSize',15,'FontName','Arial');

annotation('arrow',[2.5/9,2.5/9],[10.5/11,2.4/11],'LineWidth',1.5);
set(gca,'unit','centimeters','position',[2.5 2.5 5 8]);
set(gcf,'unit','centimeters','position',[12 8 9 11]);

