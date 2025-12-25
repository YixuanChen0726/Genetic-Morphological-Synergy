clear all;clc;close all

% Sample8
PhaCell=ElonExtreme{2,9};MusCell=ElonExtreme{5,9};IntCell=ElonExtreme{8,9};
PhaFrame=ElonExtreme{3,9};MusFrame=ElonExtreme{6,9};IntFrame=ElonExtreme{9,9};

% name dictionary
NameDic=readcell(['..\..\bin\DataSet CMap\name_dictionary.csv']);
PhaCellIndex=NameDic{find(strcmp(NameDic(:,2),PhaCell)),1};
MusCellIndex=NameDic{find(strcmp(NameDic(:,2),MusCell)),1};
IntCellIndex=NameDic{find(strcmp(NameDic(:,2),IntCell)),1};

load(['.\WT_Sample8_AllSpace.mat']);
PhaSpace=AllSpace{1,PhaFrame};MusSpace=AllSpace{1,MusFrame};IntSpace=AllSpace{1,IntFrame};

PharynxColor=[158,119,221]./255;
MuscleColor=[235,143,182]./255;
IntestineColor=[250,199,86]./255;

% Pharynx:high elongation
figure(1)
ind=find(PhaSpace==PhaCellIndex);[X,Y,Z]=ind2sub(size(PhaSpace),ind);
CellAxis=[X,Y,Z];
[SurfaceAxis,~]=boundary(CellAxis,0.95);
trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',PharynxColor,'FaceAlpha',1);hold on;
material dull;lighting gouraud;
view([-158,35]);
camlight('headlight','infinite');
axis equal off;
set(gcf,'color','k');
set(gcf, 'InvertHardCopy','off');
set(gcf,'unit','centimeters','position',[3,10,10,8]);
print(gcf,['.\Sample8_',PhaCell,'_Elongation_',num2str(PhaFrame,'%03d')],'-dpng','-r600');

% Muscle:high elongation
figure(2)
ind=find(MusSpace==MusCellIndex);[X,Y,Z]=ind2sub(size(MusSpace),ind);
CellAxis=[X,Y,Z];
[SurfaceAxis,~]=boundary(CellAxis,0.95);
trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',MuscleColor,'FaceAlpha',1);hold on;
material dull;lighting gouraud;
view([-141,25]);
camlight('headlight','infinite');
axis equal off;
set(gcf,'color','k');
set(gcf, 'InvertHardCopy','off');
set(gcf,'unit','centimeters','position',[14,10,10,8]);
print(gcf,['.\Sample8_',MusCell,'_Elongation_',num2str(MusFrame,'%03d')],'-dpng','-r600');

% Intestine:low elongation
figure(3)
ind=find(IntSpace==IntCellIndex);[X,Y,Z]=ind2sub(size(IntSpace),ind);
CellAxis=[X,Y,Z];
[SurfaceAxis,~]=boundary(CellAxis,0.95);
trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',IntestineColor,'FaceAlpha',1);hold on;
material dull;lighting gouraud;
view([-103,17]);
camlight('headlight','infinite');
axis equal off;
set(gcf,'color','k');
set(gcf, 'InvertHardCopy','off');
set(gcf,'unit','centimeters','position',[25,10,10,8]);
print(gcf,['.\Sample8_',IntCell,'_Elongation_',num2str(IntFrame,'%03d')],'-dpng','-r600');

