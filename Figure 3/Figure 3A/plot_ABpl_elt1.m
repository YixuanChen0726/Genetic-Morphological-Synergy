% plot elt-1 tree of ABpl
clear all;clc;close all;

% lineage tree
load(['..\..\bin\LineageTree.mat']);
% fate tree
load(['..\..\bin\FateClassification.mat']);
% Position Tree
load(['..\..\bin\CellPositionTree_All.mat']);
% cell gene data
load(['..\..\bin\Coasymmetry\all gene\ABpl_AllGeneSig.mat']);
ABplAll=AllGeneSigni;

% cell lifecycle
AllLifeCycle=readcell(['..\..\bin\DataSet CMap\Life Cycle\Nucleus\MeanTot_WT_Nucleus_LifeCycle.csv']);
% the last frame of 4-cell stage
ABaRow=find(strcmp(AllLifeCycle(1,:),'ABa'));ABpRow=find(strcmp(AllLifeCycle(1,:),'ABp'));
EMSRow=find(strcmp(AllLifeCycle(1,:),'EMS'));P2Row=find(strcmp(AllLifeCycle(1,:),'P2'));
StartTime=min([AllLifeCycle{3,ABaRow},AllLifeCycle{3,ABpRow},AllLifeCycle{3,EMSRow},AllLifeCycle{3,P2Row}]);
AllLifeCycle(2:end,2:end)=num2cell(cell2mat(AllLifeCycle(2:end,2:end))-StartTime);

% colorbar
c=colormap(slanCM(['BuPu']));
r1=floor(size(c,1)/3);r2=size(c,1)-r1;
c2=interp1([r2:1:size(c,1)],c(r2:size(c,1),:),[r2:0.4:size(c,1)]);
colormap([c(r1+1:r2-1,:);c2]);
BGColor=[210,210,210]./255;
LinW=1.5;MakS=1;

% ABpl and its offspring
k=find(contains(FateClas(:,1),'ABpl'));
CellList=FateClas(k,1);

% load gene data
genename='elt-1';construct='Protein';literature='5';
expdata=[genename,'_',construct,'_',literature];

% load data
load([expdata]);
% map to 0-1
RawGeneExp=[];
for CellLayer=size(LineageTree,2)-1:-1:1
    for CellIndex=1:size(LineageTree,1)
        CellName=LineageTree{CellIndex,CellLayer};
        if ~isempty(find(strcmp(CellList,CellName)))
            Exp=find(strcmp(MeanExp(:,1),CellName));
            RawGeneExp=[RawGeneExp;cell2mat(MeanExp(Exp,4))];
        end
    end
end
% maximum gene expression
[~,PS]=mapminmax(RawGeneExp',0,1);
% data list
GeneLineage=cell(size(LineageTree));
for CellLayer=size(LineageTree,2)-1:-1:1
    for CellIndex=1:size(LineageTree,1)
        CellName=LineageTree{CellIndex,CellLayer};
        q=find(strcmp(CellList,CellName));
        if ~isempty(q)
            % x position
            Xp=PositionTree{CellIndex,CellLayer};
            if all(CellName) && ~isempty(find(strcmp(AllLifeCycle(1,:),CellName)))
                Dau1Name=LineageTree{CellIndex*2-1,CellLayer+1};
                Dau2Name=LineageTree{CellIndex*2,CellLayer+1};
                Dau1Row=find(strcmp(AllLifeCycle(1,:),Dau1Name));
                Dau2Row=find(strcmp(AllLifeCycle(1,:),Dau2Name));
                % bith and death time
                MotherRow=find(strcmp(AllLifeCycle(1,:),CellName));
                BirthTime=AllLifeCycle{2,MotherRow};
                DeathTime=AllLifeCycle{3,MotherRow};
                DeathTime=max([DeathTime,AllLifeCycle{2,Dau1Row},AllLifeCycle{2,Dau2Row}]);
                % gene expression data
                Exp=find(strcmp(MeanExp(:,1),CellName));
                TimeList=cell2mat(MeanExp(Exp,2))-StartTime;
                GeneExp=mapminmax('apply',cell2mat(MeanExp(Exp,4))',PS)';
                % plot data
                % [x,birth & death, expression]
                CharaLineage=[repmat(Xp,size(TimeList)),TimeList,GeneExp;Xp,DeathTime,GeneExp(end)];
                CharaLineage=[Xp,BirthTime,CharaLineage(1,3);CharaLineage;NaN,NaN,NaN];
                % save format
                GeneLineage{CellIndex,CellLayer}=CharaLineage;
            end
        end
    end
end

% background
for CellLayer=size(GeneLineage,2):-1:1
    for CellIndex=1:1:size(GeneLineage,1)
        Data=GeneLineage{CellIndex,CellLayer};
        if ~isempty(Data)
            Xp=Data(1,1);
            BirthTime=Data(1,2);
            DeathTime=Data(end-1,2);
            plot([Xp,Xp],[BirthTime,DeathTime],'Color',BGColor,'LineWidth',LinW);hold on;
            MData=GeneLineage{ceil(CellIndex/2),CellLayer-1};
            if ~isempty(MData)
                MXp=MData(1,1);
                % connect line
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
set(get(c1,'ylabel'),'string',['\it',genename,'\rm'],'fontsize',15);
clim([0,1]);set(c1,'ticks',0:0.2:1);

% axis range
xrange=[];
for CellIndex=1:size(GeneLineage,1)
    for CellLayer=1:size(GeneLineage,2)
        Data=GeneLineage{CellIndex,CellLayer};
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

% color bar
function colorList=slanCM(type,num)
if nargin<2
    num=256;
end
if nargin<1
    type='';
end

slanCM_Data=load('slanCM_Data.mat');
CList_Data=[slanCM_Data.slandarerCM(:).Colors];

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
