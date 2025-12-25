% Sample7,Daa,107
clear all;clc;close all;

% Name Dictionary
NameDic=readcell(['..\..\bin\name_dictionary.csv']);
CellIndex=NameDic{find(strcmp(NameDic(:,2),'Daa')),1};
% Lineage Tree
load(['..\..\bin\LineageTree.mat']);
% Colorbar
load(['..\..\bin\ColorTree2.mat']);
[a,b]=find(strcmp(LineageTree,'Daa'));
Color=ColorTree{a,b};
Color1=[242,192,215]./255;
Color2=[181,59,114]./255;

% plot
AxisRotation=[1,-1,-1];
Space=niftiread(['.\WT_Sample7_107_segCell.nii.gz']);
ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
CellAxis=[AxisRotation(1,1)*X,AxisRotation(1,2)*Y,AxisRotation(1,3)*Z];
[SurfaceAxis,~]=boundary(CellAxis,0.95);

% [a,a axis,b,b axis,c, caxis]
load(['..\..\Morphology Descriptor\OBB box\WT_Sample7_OBBbox.mat']);
IndexRow=find(cell2mat(OBBbox(:,3))==CellIndex);TimeRow=find(cell2mat(OBBbox(:,4))==107);
OBBRow=intersect(IndexRow,TimeRow);
a=OBBbox{OBBRow,5};b=OBBbox{OBBRow,7};c=OBBbox{OBBRow,9};
aspace=OBBbox{OBBRow,6};bspace=OBBbox{OBBRow,8};cspace=OBBbox{OBBRow,10};

% centroid
CellAxis=[CellAxis(:,1)-mean(CellAxis(:,1)),CellAxis(:,2)-mean(CellAxis(:,2)),CellAxis(:,3)-mean(CellAxis(:,3))];
minX=min(CellAxis(:,1));maxX=max(CellAxis(:,1));
minY=min(CellAxis(:,2));maxY=max(CellAxis(:,2));
minZ=min(CellAxis(:,3));maxZ=max(CellAxis(:,3));
[SurfaceAxis,VolumeTemp]=boundary(CellAxis,0.98);

figure(1)
% start:[minX,minY,minZ]
aAxis=aspace-repmat((aspace(1,:)-[minX+10,minY,maxZ+10]),2,1);
bAxis=bspace-repmat((bspace(2,:)-[minX+10,minY,maxZ+10]),2,1);
cAxis=cspace-repmat((cspace(1,:)-[minX+10,minY,maxZ+10]),2,1);

trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',Color1,'FaceAlpha',1);hold on;
plot3(aAxis(:,1),aAxis(:,2),aAxis(:,3),'Color','k','LineWidth',2);hold on;
plot3(bAxis(:,1),bAxis(:,2),bAxis(:,3),'Color','k','LineWidth',2);hold on;
plot3(cAxis(:,1),cAxis(:,2),cAxis(:,3),'Color','k','LineWidth',2);hold on;
text(aAxis(2,1)+2,aAxis(2,2)-5,aAxis(2,3)-5,'\ita','Color','k','FontSize',20,'Fontname','Arial','FontWeight','bold');
text(bAxis(1,1)+4,bAxis(1,2)+5,bAxis(1,3)+5,'\itb','Color','k','FontSize',20,'Fontname','Arial','FontWeight','bold');
text(cAxis(2,1)+2,cAxis(2,2)-2,cAxis(2,3)-2,'\itc','Color','k','FontSize',20,'Fontname','Arial','FontWeight','bold');

axis off;box off;axis equal;
view(-22,-21);
camlight('headlight','infinite');
lighting gouraud;
material([0.6,0.4,0,1]);
set(gcf,'color','w');
set(gcf, 'InvertHardCopy','off');
set(gcf,'unit','centimeters','position',[10,8,15,10]);
set(gca,'unit','centimeters','position',[1,1,13,8]);

