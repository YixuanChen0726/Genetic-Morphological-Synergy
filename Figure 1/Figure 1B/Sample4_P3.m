% P Lineage:asymmetric division
clear all;clc;close all;

Mothername='P3';
dau1name='D';
dau2name='P4';
load('..\..\bin\LineageTree.mat');
load('..\..\bin\ColorTree2.mat');
AxisRotation=[1,-1,-1;1,1,1;1,1,1;-1,-1,1;1,-1,-1;-1,1,-1;1,-1,-1;-1,-1,1];

% color list
[MRow,MColumn]=find(strcmp(LineageTree,Mothername));
[D1Row,D1Column]=find(strcmp(LineageTree,dau1name));
[D2Row,D2Column]=find(strcmp(LineageTree,dau2name));
ColorMother=ColorTree{MRow,MColumn};
Colordau1=ColorTree{D1Row,D1Column};
Colordau2=ColorTree{D2Row,D2Column};

% load name dictionary
NameDic=readcell(['..\..\bin\name_dictionary.csv']);
CellIndexMother=find(strcmp(NameDic(:,2),Mothername))-1;
CellIndexdau1=find(strcmp(NameDic(:,2),dau1name))-1;
CellIndexdau2=find(strcmp(NameDic(:,2),dau2name))-1;

SampleNum=4;
% P3 lifecycle
LifeCycle=readcell(['..\..\Morphology Descriptor\LifeCycle\WT_Sample',num2str(SampleNum),'_Membrane_LifeCycle.csv']);
NameIndexMother=find(strcmp(LifeCycle(1,:),Mothername));
DeathTimeMother=LifeCycle{4,NameIndexMother};
NameIndexdau1=find(strcmp(LifeCycle(1,:),dau1name));
BirthTimedau1=LifeCycle{2,NameIndexdau1};
NameIndexdau2=find(strcmp(LifeCycle(1,:),dau2name));
BirthTimedau2=LifeCycle{2,NameIndexdau2};
BirthTimedau=max([BirthTimedau1,BirthTimedau2]);

% load nii.gz Segmemtation Data
load(['.\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
% Mother's space
Space=AllSpace{1,DeathTimeMother};
ind=find(Space==CellIndexMother);[X,Y,Z]=ind2sub(size(Space),ind);
CellAxisMother=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
% Daughter' space
Space=AllSpace{1,BirthTimedau};
ind=find(Space==CellIndexdau1);[X,Y,Z]=ind2sub(size(Space),ind);
CellAxisdau1=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
ind=find(Space==CellIndexdau2);[X,Y,Z]=ind2sub(size(Space),ind);
CellAxisdau2=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
clear AllSpace

angle=[-190,10];
t = tiledlayout(12,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

figure(1)
CellAxisMother=[CellAxisMother(:,1)-mean(CellAxisMother(:,1)),CellAxisMother(:,2)-mean(CellAxisMother(:,2)),CellAxisMother(:,3)-mean(CellAxisMother(:,3))];
[SurfaceAxis,VolumeTemp]=boundary(CellAxisMother,0.96);
trisurf(SurfaceAxis,CellAxisMother(:,1),CellAxisMother(:,2),CellAxisMother(:,3),'EdgeColor','none','FaceColor',ColorMother,'FaceAlpha',1);hold on;
view(angle(1,:));
camlight('headlight','infinite');
material dull;lighting phong;
axis off;box off;axis square;axis equal;
axis([min(CellAxisMother(:,1))-2,max(CellAxisMother(:,1))+2,min(CellAxisMother(:,2))-2,max(CellAxisMother(:,2))+2,min(CellAxisMother(:,3))-2,max(CellAxisMother(:,3))+2]);
set(gcf,'color','k');
set(gcf, 'InvertHardCopy','off');
set(gca,'unit','centimeters','position',[1 1 13 8]);
set(gcf,'unit','centimeters','position',[5 5 15 10]);

figure(2)
CellAxis=[CellAxisdau1;CellAxisdau2];
CellAxisdau1=[CellAxisdau1(:,1)-mean(CellAxis(:,1)),CellAxisdau1(:,2)-mean(CellAxis(:,2)),CellAxisdau1(:,3)-mean(CellAxis(:,3))];
CellAxisdau2=[CellAxisdau2(:,1)-mean(CellAxis(:,1)),CellAxisdau2(:,2)-mean(CellAxis(:,2)),CellAxisdau2(:,3)-mean(CellAxis(:,3))];
CellAxis=[CellAxis(:,1)-mean(CellAxis(:,1)),CellAxis(:,2)-mean(CellAxis(:,2)),CellAxis(:,3)-mean(CellAxis(:,3))];
[SurfaceAxisdau1,~]=boundary(CellAxisdau1,0.96);
trisurf(SurfaceAxisdau1,CellAxisdau1(:,1),CellAxisdau1(:,2),CellAxisdau1(:,3),'EdgeColor','none','FaceColor',Colordau1,'FaceAlpha',1);hold on;
[SurfaceAxisdau2,~]=boundary(CellAxisdau2,0.96);
trisurf(SurfaceAxisdau2,CellAxisdau2(:,1),CellAxisdau2(:,2),CellAxisdau2(:,3),'EdgeColor','none','FaceColor',Colordau2,'FaceAlpha',1);hold on;
view(angle(1,:));
camlight('headlight','infinite');
material dull;lighting phong;
axis off;box off;axis square;axis equal;
axis([min(CellAxis(:,1))-2,max(CellAxis(:,1))+2,min(CellAxis(:,2))-2,max(CellAxis(:,2))+2,min(CellAxis(:,3))-2,max(CellAxis(:,3))+2]);

set(gcf,'color','k');
set(gcf, 'InvertHardCopy','off');
set(gca,'unit','centimeters','position',[1 1 13 8]);
set(gcf,'unit','centimeters','position',[22 5 15 10]);
