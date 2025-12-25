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

% 画图
AxisRotation=[1,-1,-1];
Space=niftiread(['.\WT_Sample7_107_segCell.nii.gz']);
ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
CellAxis=[AxisRotation(1,1)*X,AxisRotation(1,2)*Y,AxisRotation(1,3)*Z];
[SurfaceAxis,~]=boundary(CellAxis,0.95);

figure(1)
h=trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor',Color2,'LineWidth',0.2,'FaceColor',Color1,'FaceAlpha',1);hold on;
axis off;box off;axis equal;
view(-39,-34);
camlight('headlight','infinite');
lighting gouraud;
material([0.6,0.4,0,1]);
set(gca,'color','none');
set(gcf,'color','w');
set(gcf, 'InvertHardCopy','off');
set(gcf,'unit','centimeters','position',[10,8,15,10]);
set(gca,'unit','centimeters','position',[1,1,13,8]);
