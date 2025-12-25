% the developmental images of the embryo, colored by elt-1 expression level
clear all;clc;close all;

% the rotation of 8 embryos
AxisRotation=[1,-1,-1;1,1,1;1,1,1;-1,-1,1;1,-1,-1;-1,1,-1;1,-1,-1;-1,-1,1];
viewangle=[-90,90;90,-90;-90,180;-90,0;0,180;180,-180];
% name dictionary
NameDic=readcell(['..\..\bin\DataSet CMap\name_dictionary.csv']);
% lineange tree
load('..\..\bin\LineageTree.mat');
% fate destribution
load('..\..\binFateClassification.mat');

% colorbar
c=colormap(slanCM('RdBu'));
% interp
r1=floor(size(c,1)/3);r2=size(c,1)-ceil(size(c,1)/3);
r3=floor(size(c,1)/7);r4=size(c,1)-ceil(size(c,1)/7);
midr=(c(r1,:)+c(r2,:))./2;
c0=interp1([1:1:3],[c(r1,:);midr;c(r2,:)],[1:0.05:3]);
c3=interp1([1:1:r3],[c(1:r3,:)],[1:0.16:r3]);
c4=interp1([r4:1:size(c,1)],[c(r4:end,:)],[r4:0.6:size(c,1)]);
newc=[c3;c(r3+1:r1-1,:);c0;c(r2+1:r4-1,:);c4];
p=colormap(flipud(newc));
colormap(p);

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

% shape descriptor data
load(['..\..\bin\DataSet CMap\Elongation Ratio\WT_Sample',num2str(SampleNum),'_Elongation.mat']);
shapetable=Elongation(:,[2,4,5]);
RawShape=cell2mat(shapetable(:,3));
% mapr to colorbar
[G,PS]=mapminmax(RawShape',1,size(p,1));

% plot
for Frame=1:length(AllSpace)
    Space=AllSpace{1,Frame};
    ExistCell=double(unique(Space));
    ExistCell(find(ExistCell==0))=[];
    CellNumber=length(ExistCell);
    % plot frame
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
        Row1=find(strcmp(shapetable(:,1),CellName));
        Row2=find(cell2mat(shapetable(:,2))==Frame);
        Row=intersect(Row1,Row2);
        if ~isempty(Row)
            Exp=shapetable{Row,3};
            ShapeExp=mapminmax('apply',Exp,PS);
            Color=p(round(ShapeExp),:);
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
    % background parameters
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
