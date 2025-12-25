% mesoderm: diameter sphericity
clear all; clc; close all;

% Name Dictionary
NameDic=readcell(['..\..\..\bin\DataSet CMap\name_dictionary.csv']);
NeuronColor=[126,186,210]./255;
PharynxColor=[158,119,221]./255;
MuscleColor=[235,143,182]./255;
IntestineColor=[250,199,86]./255;
SkinColor=[147,202,80]./255;
MesoColor=[216,108,88]./255;
AxisRotation=[1,-1,-1;1,1,1;1,1,1;-1,-1,1;1,-1,-1;-1,1,-1;1,-1,-1;-1,-1,1];

load(['..\Mesoderm_ExtremeShapeList.mat']);
% diameter sphericity
HighData=ExtremeFateList{2,2};
LowData=ExtremeFateList{3,2};

viewangle=[80,13;119,12;13,9;103,-10;-90,-60;
    196,35;13,7;63,7;40,3;68,5];
for Num=1:10
    figure(Num)
    DSphHigh=HighData(Num,:);
    SampleNum=DSphHigh{1,1};CellName=DSphHigh{1,2};Frame=DSphHigh{1,3};
    CellIndex=NameDic{find(strcmp(NameDic(:,2),CellName)),1};
    clear AllSpace;
    load(['..\..\..\bin\CMap Data\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
    Space=AllSpace{1,Frame};
    ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
    CellAxis=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
    [SurfaceAxis,~]=boundary(CellAxis,0.95);
    h=trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',MesoColor,'FaceAlpha',1);hold on;
    axis off;box off;axis equal;
    view(viewangle(Num,:));
    camlight('headlight','infinite');
    material dull;lighting gouraud;
    title(CellName,'FontSize',15,'FontName','Arial','Color','w');
    set(gcf,'color','k');
    set(gcf, 'InvertHardCopy','off');
    set(gcf,'unit','centimeters','position',[10,5,8,6]);
    set(gca,'unit','centimeters','position',[2,1,5,3.5]);
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_DSphHigh',num2str(Num)],'-dpng','-r600');
end

viewangle=[223,-2;-114,-8;-118,6;-163,-3;16,34;
    -96,18;19,8;9,5;-72,12;-74,6];
for Num=1:10
    figure(Num+10)
    DSphLow=LowData(Num,:);
    SampleNum=DSphLow{1,1};CellName=DSphLow{1,2};Frame=DSphLow{1,3};
    CellIndex=NameDic{find(strcmp(NameDic(:,2),CellName)),1};
    clear AllSpace;
    load(['..\..\..\bin\CMap Data\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
    Space=AllSpace{1,Frame};
    ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
    CellAxis=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
    [SurfaceAxis,~]=boundary(CellAxis,0.95);
    h=trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',MesoColor,'FaceAlpha',1);hold on;
    axis off;box off;axis equal;
    view(viewangle(Num,:));
    camlight('headlight','infinite');
    material dull;lighting gouraud;
    title(CellName,'FontSize',15,'FontName','Arial','Color','w');
    set(gcf,'color','k');
    set(gcf, 'InvertHardCopy','off');
    set(gcf,'unit','centimeters','position',[10,5,8,6]);
    set(gca,'unit','centimeters','position',[2,1,5,3.5]);
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_DSphLow',num2str(Num)],'-dpng','-r600');
end

