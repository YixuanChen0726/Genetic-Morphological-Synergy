% mesoderm: general sphericity
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
% general sphericity
HighData=ExtremeFateList{2,1};
LowData=ExtremeFateList{3,1};

viewangle=[95,5;-196,10;107,21;86,-17;-100,-84;
    192,77;86,5;-104,5;108,16;170,15];
for Num=1:10
    figure(Num)
    GSphHigh=HighData(Num,:);
    SampleNum=GSphHigh{1,1};CellName=GSphHigh{1,2};Frame=GSphHigh{1,3};
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
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_GSphHigh',num2str(Num)],'-dpng','-r600');
end

viewangle=[180,-2;119,-14;-168,-78;-52,-9;58,-13;
    25,4;-100,77;187,-44;-60,-56;-53,2];
for Num=1:10
    figure(Num+10)
    GSphLow=LowData(Num,:);
    SampleNum=GSphLow{1,1};CellName=GSphLow{1,2};Frame=GSphLow{1,3};
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
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_GSphLow',num2str(Num)],'-dpng','-r600');
end

