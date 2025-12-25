% plot lineage tree colored by Wilson Flatness Index
clear all;clc;close all;

% lineage tree
load(['..\..\bin\LineageTree.mat']);
% fate tree
load(['..\..\bin\FateClassification.mat']);
% Position Tree
load(['..\..\bin\CellPositionTree_All.mat']);
% cell shape data
load(['..\..\bin\DataSet CMap\Shape Lineage Tree\WilsonFlatnessIndex_AllLineage.mat']);

% cell life cycle
AllLifeCycle=readcell(['..\..\bin\DataSet CMap\Life Cycle\Nucleus\MeanTot_WT_Nucleus_LifeCycle.csv']);
% the last frame of 4-cell stage
ABaRow=find(strcmp(AllLifeCycle(1,:),'ABa'));ABpRow=find(strcmp(AllLifeCycle(1,:),'ABp'));
EMSRow=find(strcmp(AllLifeCycle(1,:),'EMS'));P2Row=find(strcmp(AllLifeCycle(1,:),'P2'));
StartTime=min([AllLifeCycle{3,ABaRow},AllLifeCycle{3,ABpRow},AllLifeCycle{3,EMSRow},AllLifeCycle{3,P2Row}]);
AllLifeCycle(2:end,2:end)=num2cell(cell2mat(AllLifeCycle(2:end,2:end))-StartTime);

% Capa: ventral muscle
CapaVen={'Capaa';'Capaaa';'Capaap';'Capapp';'Capaaaa';'Capaaap';'Capaapa';'Capaapp';'Capappp';'Capappa'};
% Capp: dorsal muscle
CappDor={'Cappa';'Cappaa';'Capppa';'Cappap';'Cappaaa';'Cappaap';'Capppaa';'Capppap';'Cappapa';'Cappapp';'Cappppd'};
PartialCellList=[CapaVen;CappDor];

% color
c=colormap(slanCM('RdBu'));
c=flipud(c);
% interp
r1=floor(size(c,1)/3.5);r2=size(c,1)-r1;
c1=interp1([1:1:r1],c(1:r1,:),[1:0.38:r1]);
c2=interp1([r2:1:size(c,1)],c(r2:size(c,1),:),[r2:0.2:size(c,1)]);
p=colormap([c1;c(r1+1:r2-1,:);c2]);

k=find(contains(FateClas(:,1),'Cap'));
CellList=FateClas(k,1);
CellList=intersect(CellList,AllLifeCycle(1,2:end));

ShapeLineage=cell(size(LineageData));
for i=1:length(CellList)
    CellName=CellList{i,1};
    [a,b]=find(strcmp(LineageTree,CellName));
    % cell position
    Position=PositionTree{a,b};
    CellData=LineageData{a,b};
    % modify life cycle
    MotherRow=find(strcmp(AllLifeCycle(1,:),CellName));
    if ~isempty(MotherRow)
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
    else
        GrandName=LineageTree{ceil(a/2),b-1};
        GrandRow=find(strcmp(AllLifeCycle(1,:),GrandName));
        BirthTime=AllLifeCycle{3,GrandRow};
        DeathTime=max(cell2mat(AllLifeCycle(3,2:end)));
    end
    % map cell life
    if ~isempty(CellData)
        NewLife=mapminmax(CellData(1:end-1,2)',BirthTime,DeathTime);
        % [x, time point, mean data, std data]
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
        CellName=LineageTree{CellIndex,CellLayer};
        if ~isempty(Data)
            Xp=Data(1,1);
            BirthTime=Data(1,2);
            DeathTime=Data(end-1,2);
            appear=find(strcmp(PartialCellList,CellName));
            plot([Xp,Xp],[BirthTime,DeathTime],'Color',BGColor,'LineWidth',LinW);hold on;
            MData=ShapeLineage{ceil(CellIndex/2),CellLayer-1};
            if ~isempty(MData)
                MXp=MData(1,1);
                plot([Xp,MXp],[BirthTime,BirthTime],'Color',BGColor,'LineWidth',LinW);hold on;
            end
            if ~isnan(Data(1,3)) && ~isempty(appear)
                patch(Data(:,1),Data(:,2),Data(:,3),'EdgeColor','interp','Marker','none','MarkerSize',MakS,'MarkerFaceColor','flat','LineWidth',LinW);hold on;
            end
            if ~isempty(MData) && ~isempty(appear)
                MXp=MData(1,1);
                patch([Xp,MXp],[Data(1,2),Data(1,2)],[Data(1,3),Data(1,3)],'EdgeColor','interp','Marker','none','MarkerSize',MakS,'MarkerFaceColor','flat','LineWidth',LinW);hold on;
            end
        end
    end
end

c1=colorbar;
set(get(c1,'ylabel'),'string','Wilson Flatness Index','fontsize',15);
set(c1,'ticks',0.05:0.1:2);

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

% mark Lineage
[Row,Column]=find(strcmp(LineageTree,'Cap'));
X=ShapeLineage{Row,Column}(1,1);
Y=ShapeLineage{Row,Column}(1,2);
text(X,Y-18,'Cap','FontSize',15,'FontName','Arial','HorizontalAlignment','center');

axis([min(xrange)-2,max(xrange)+2,30,320]);
yticks(0:40:280);box off;
ylabel({'\it t\rm (min)'},'Rotation',0,'Position',[min(xrange)-7 330]);
set(gca,'YDir','reverse','YAxisLocation','left','ycolor','k');
set(gca,'xtick',[],'XAxisLocation','top','xcolor','none');
set(gca,'LineWidth',1,'FontSize',15,'FontName','Arial');

% partial mean level/asymmetry/p value
load(['..\Figure 6C\Cap_partial_shapesig.mat']);
n=find(strcmp(ShapeSigni(1,:),'Wilson Flatness index'));
meandata=ShapeSigni{4,n};
pvalue=ShapeSigni{5,n};pindex=ceil(-log10(pvalue));
asym=ShapeSigni{6,n};

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

plot([Dau1range(1,1),Dau1range(2,1)],[300,300],'-','LineWidth',LinW,'Color','k');hold on;
plot([Dau2range(1,1),Dau2range(2,1)],[300,300],'-','LineWidth',LinW,'Color','k');hold on;

if meandata(1,1)>meandata(2,1)
    text(mean(Dau1range),315,'high','FontSize',15,'FontName','Arial','HorizontalAlignment','center');
    text(mean(Dau2range),315,'low','FontSize',15,'FontName','Arial','HorizontalAlignment','center');
else
    text(mean(Dau1range),315,'low','FontSize',15,'FontName','Arial','HorizontalAlignment','center');
    text(mean(Dau2range),315,'high','FontSize',15,'FontName','Arial','HorizontalAlignment','center');
end
text(mean([Dau1range(1,1);Dau2range(2,1)]),345,['\itp\rm = ',num2str(pvalue*10^pindex,'%.2f'),'\times10^{-',num2str(pindex),'}'],'FontSize',15,'FontName','Arial','HorizontalAlignment','center');

axis([min(xrange)-2,max(xrange)+2,30,320]);
yticks(0:40:280);box off;
ylabel({'\it t\rm (min)'},'Rotation',0,'Position',[min(xrange)-7 330]);
set(gca,'YDir','reverse','YAxisLocation','left','ycolor','k');
set(gca,'xtick',[],'XAxisLocation','top','xcolor','none');
set(gca,'LineWidth',1,'FontSize',15,'FontName','Arial');

annotation('arrow',[2.5/10.5,2.5/10.5],[10/11,2/11],'LineWidth',1.5);
set(gca,'unit','centimeters','position',[2.5 2 5 8]);
set(gcf,'unit','centimeters','position',[12 8 10.5 11]);

% save image
saveas(gcf,['Cap_partial_WilsonFlatnessIndex.svg'],'svg');
print(gcf,['Cap_partial_WilsonFlatnessIndex'],'-dpng','-r600');

% colorbar function
function colorList=slanCM(type,num)
if nargin<2
    num=256;
end
if nargin<1
    type='';
end

slanCM_Data=load('..\..\Figure 2\slanCM_Data.mat');
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
