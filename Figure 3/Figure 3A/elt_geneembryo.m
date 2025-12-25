% the developmental images of the embryo, colored by elt-1 expression level
clear all;clc;close all;

% the rotation of 8 embryos
AxisRotation=[1,-1,-1;1,1,1;1,1,1;-1,-1,1;1,-1,-1;-1,1,-1;1,-1,-1;-1,-1,1];
viewangle=[-90,90;90,-90;-90,180;-90,0;0,180;180,-180];
% name dictionary
NameDic=readcell(['..\..\bin\DataSet CMap\name_dictionary.csv']);
% lineage tree
load('..\..\bin\LineageTree.mat');
% fate destribution
load('..\..\bin\FateClassification.mat');

% load gene edata
genename='elt-1';construct='Protein';literature='5';
% WorkSpace_5_151_2_1.csv

c=colormap(slanCM(['BuPu']));
p=colormap(c([floor(size(c,1)/3):size(c,1)],:));

OtherColor=[200,200,200]./255;

figure
SampleNum=7;
% data of all the membane
clear AllSpace;
load(['.\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
mkdir(['.\Sample',num2str(SampleNum)]);
FrameNum=length(AllSpace);
% axis range
meanAxis=[];maxAxis=[];minAxis=[];
for Frame=1:length(AllSpace)
    Space=AllSpace{1,Frame};
    ExistCell=unique(Space);
    ExistCell(find(ExistCell==0))=[];
    AllCellAxis=[];
    for Vari=1:length(ExistCell)
        CellIndex=ExistCell(Vari);
        ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
        AllCellAxis=[AllCellAxis;AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
    end
    % mean/max/min data
    meanAxis=[meanAxis;mean(AllCellAxis(:,1)),mean(AllCellAxis(:,2)),mean(AllCellAxis(:,3))];
    maxAxis=[maxAxis;max(AllCellAxis(:,1)),max(AllCellAxis(:,2)),max(AllCellAxis(:,3))];
    minAxis=[minAxis;min(AllCellAxis(:,1)),min(AllCellAxis(:,2)),min(AllCellAxis(:,3))];
end
% centroid
centroid=mean(meanAxis,1);
CmaxAxis=max(maxAxis,[],1)-centroid;CminAxis=min(minAxis,[],1)-centroid;
axisrange=[CminAxis;CmaxAxis]';
clear AllCellAxis;

% gene data:WorkSpace_5_102_2_1.csv
genetable=readcell(['.\WorkSpace_5_102_2_',num2str(SampleNum),'.csv']);
genetable=genetable(2:end,[1,2,4]);
RawGeneExp=cell2mat(genetable(:,3));
% map to 0-1
[G,PS]=mapminmax(RawGeneExp',1,size(p,1));

% plot
for Frame=1:length(AllSpace)
    Space=AllSpace{1,Frame};
    ExistCell=double(unique(Space));
    ExistCell(find(ExistCell==0))=[];
    CellNumber=length(ExistCell);
    clf;
    set(gcf,'unit','centimeters','position',[2 8 40 12]);
    t=tiledlayout(1,6);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    titlestr={'A-P/D-V(L)','A-P/D-V(R)','A-P/L-R(D)','A-P/L-R(V)','D-V/L-R(A)','D-V/L-R(P)'};
    % diagram
    ax=nexttile(1);
    % plot by cells
    for Vari=1:length(ExistCell)
        CellIndex=ExistCell(Vari,1);
        % cell index
        CellName=NameDic{CellIndex+1,2};
        alpha=1;Color=OtherColor;
        Row1=find(strcmp(genetable(:,1),CellName));
        Row2=find(cell2mat(genetable(:,2))==Frame);
        Row=intersect(Row1,Row2);
        if ~isempty(Row)
            Exp=genetable{Row,3};
            if ~isnan(Exp)
                GeneExp=mapminmax('apply',Exp,PS);
                Color=p(round(GeneExp),:);
            end
        end
        ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
        CellAxis=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
        CellAxis=CellAxis-centroid;
        [SurfaceAxis,~]=boundary(CellAxis,0.95);
        h=trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',Color,'FaceAlpha',alpha);hold on;
    end
    axis off;box off;axis equal;
    axis([axisrange(1,1)-20,axisrange(1,2)+20,axisrange(2,1)-20,axisrange(2,2)+20,axisrange(3,1)-20,axisrange(3,2)+20]);
    % copy image
    for i=2:6
        axc=copyobj(ax,t);
        axc(1).Layout.Tile=i;
    end
    % modify image parameters
    for i=1:6
        nexttile(i)
        material dull;lighting gouraud;
        view(viewangle(i,:));
        camlight('headlight','infinite');
    end
    title(t,{['Cell Number = ',num2str(CellNumber)],[],['  ',titlestr{1,1},'                ',titlestr{1,2},'                ',titlestr{1,3},'                ',titlestr{1,4},'                ',titlestr{1,5},'                ',titlestr{1,6},'  ']},'FontSize',18,'FontName','Arial','Color','w');
    % background parameter
    set(gcf,'color','k');
    set(gcf, 'InvertHardCopy','off');
    set(t,'unit','centimeters','position',[2,1,37,7]);
    print(gcf,['.\Sample',num2str(SampleNum),'\Sample',num2str(SampleNum),'_',num2str(Frame,'%03d')],'-dpng','-r600');
    clf;
end

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
