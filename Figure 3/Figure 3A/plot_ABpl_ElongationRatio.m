% plot Elongation Ratio tree of ABpl
clear all;clc;close all;

% lineage tree
load(['..\..\bin\LineageTree.mat']);
% fate distribution
load(['..\..\bin\FateClassification.mat']);
% Position Tree
load(['..\..\bin\CellPositionTree_All.mat']);
% cell shape data
load(['..\..\bin\DataSet CMap\Shape Lineage Tree\ElongationRatio_AllLineage.mat']);

% cell nucleus life
AllLifeCycle=readcell(['..\..\bin\DataSet CMap\Life Cycle\Nucleus\MeanTot_WT_Nucleus_LifeCycle.csv']);
% the last frame of 4-cell stage
ABaRow=find(strcmp(AllLifeCycle(1,:),'ABa'));ABpRow=find(strcmp(AllLifeCycle(1,:),'ABp'));
EMSRow=find(strcmp(AllLifeCycle(1,:),'EMS'));P2Row=find(strcmp(AllLifeCycle(1,:),'P2'));
StartTime=min([AllLifeCycle{3,ABaRow},AllLifeCycle{3,ABpRow},AllLifeCycle{3,EMSRow},AllLifeCycle{3,P2Row}]);
AllLifeCycle(2:end,2:end)=num2cell(cell2mat(AllLifeCycle(2:end,2:end))-StartTime);

% colorbar
c=colormap(slanCM('RdBu'));
% interp
r1=floor(size(c,1)/3);r2=size(c,1)-ceil(size(c,1)/3);
r3=floor(size(c,1)/7);r4=size(c,1)-ceil(size(c,1)/7);
midr=(c(r1,:)+c(r2,:))./2;
c0=interp1([1:1:3],[c(r1,:);midr;c(r2,:)],[1:0.05:3]);
c3=interp1([1:1:r3],[c(1:r3,:)],[1:0.1:r3]);
c4=interp1([r4:1:size(c,1)],[c(r4:end,:)],[r4:0.45:size(c,1)]);
newc=[c3;c(r3+1:r1-1,:);c0;c(r2+1:r4-1,:);c4];
p=colormap(newc);
colormap(flipud(p));

% ABpl sublineage
k=find(contains(FateClas(:,1),'ABpl'));
CellList=FateClas(k,1);
CellList=intersect(CellList,AllLifeCycle(1,2:end));

ShapeLineage=cell(size(LineageData));
for i=1:length(CellList)
    CellName=CellList{i,1};
    [a,b]=find(strcmp(LineageTree,CellName));
    % cell position
    Position=PositionTree{a,b};
    CellData=LineageData{a,b};
    % modify cell lifecycle
    MotherRow=find(strcmp(AllLifeCycle(1,:),CellName));
    if ~isempty(MotherRow)
        BirthTime=AllLifeCycle{2,MotherRow};
        % daughter
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
    else
        GrandName=LineageTree{ceil(a/2),b-1};
        GrandRow=find(strcmp(AllLifeCycle(1,:),GrandName));
        BirthTime=AllLifeCycle{3,GrandRow};
        DeathTime=max(cell2mat(AllLifeCycle(3,2:end)));
    end
    % mapping
    if ~isempty(CellData)
        NewLife=mapminmax(CellData(1:end-1,2)',BirthTime,DeathTime);
        % [x, timepoint, meanData, stdData]
        CharaLineage=[CellData(:,1),[NewLife';NaN],CellData(:,3)];
    else
        CharaLineage=[repmat(Position,3,1),[BirthTime;DeathTime;NaN],nan(3,1)];
    end
    ShapeLineage{a,b}=CharaLineage;
end

% background
BGColor=[210,210,210]./255;
LinW=1.5;MakS=1;
for CellLayer=size(ShapeLineage,2):-1:1
    for CellIndex=1:1:size(ShapeLineage,1)
        Data=ShapeLineage{CellIndex,CellLayer};
        if ~isempty(Data)
            Xp=Data(1,1);
            BirthTime=Data(1,2);
            DeathTime=Data(end-1,2);
            plot([Xp,Xp],[BirthTime,DeathTime],'Color',BGColor,'LineWidth',LinW);hold on;
            % mother exist or not
            MData=ShapeLineage{ceil(CellIndex/2),CellLayer-1};
            if ~isempty(MData)
                MXp=MData(1,1);
                % connect mother and daughter
                plot([Xp,MXp],[BirthTime,BirthTime],'Color',BGColor,'LineWidth',LinW);hold on;
            end
            if ~isnan(Data(1,3))
                patch(Data(:,1),Data(:,2),Data(:,3),'EdgeColor','interp','Marker','none','MarkerSize',MakS,'MarkerFaceColor','flat','LineWidth',LinW);hold on;
            end
            if ~isempty(MData)
                MXp=MData(1,1);
                patch([Xp,MXp],[Data(1,2),Data(1,2)],[Data(1,3),Data(1,3)],'EdgeColor','interp','Marker','none','MarkerSize',MakS,'MarkerFaceColor','flat','LineWidth',LinW);hold on;
            end
        end
    end
end

c1=colorbar;
set(get(c1,'ylabel'),'string','Elongation Ratio','fontsize',15);
set(c1,'ticks',0:0.2:3);

% axis range
xrange=[];
for CellIndex=1:size(ShapeLineage,1)
    for CellLayer=1:size(ShapeLineage,2)
        Data=ShapeLineage{CellIndex,CellLayer};
        if ~isempty(Data)
            Position=Data(1,1);
            xrange=[xrange;Position];
        end
    end
end

axis([min(xrange)-17,max(xrange)+17,-5,310]);
yticks(0:40:280);box off;
ylabel({'\it t\rm (min)'},'Rotation',0,'Position',[min(xrange)-55 310]);
set(gca,'YDir','reverse','YAxisLocation','left','ycolor','k');
set(gca,'xtick',[],'XAxisLocation','top','xcolor','none');
set(gca,'LineWidth',1,'FontSize',15,'FontName','Arial');

annotation('arrow',[2.2/21,2.2/21],[19/20,1/20],'LineWidth',1.5);
set(gca,'unit','centimeters','position',[2.2 1 16 18]);
set(gcf,'unit','centimeters','position',[12 4 21 20]);

% color function
function colorList=slanCM(type,num)
if nargin<2
    num=256;
end
if nargin<1
    type='';
end

slanCM_Data=load('.\slanCM_Data.mat');
CList_Data=[slanCM_Data.slandarerCM(:).Colors];
disp(slanCM_Data.author);

if isnumeric(type)
    Cmap=CList_Data{type};
else
    Cpos=strcmpi(type,slanCM_Data.fullNames);
    Cmap=CList_Data{Cpos};
end

Ci=1:256;Cq=linspace(1,256,num);
colorList=[interp1(Ci,Cmap(:,1),Cq,'linear')',...
    interp1(Ci,Cmap(:,2),Cq,'linear')',...
    interp1(Ci,Cmap(:,3),Cq,'linear')'];
end
