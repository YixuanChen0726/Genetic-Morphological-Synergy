% Mesoderm: spreading index
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
% Spreading index
HighData=ExtremeFateList{2,6};
LowData=ExtremeFateList{3,6};

viewangle=[66,12;-177,31;85,-17;133,6;-168,-16;
    139,18;0,64;34,22;123,-9;177,23];
for Num=1:10
    figure(Num)
    SpreadHigh=HighData(Num,:);
    SampleNum=SpreadHigh{1,1};CellName=SpreadHigh{1,2};Frame=SpreadHigh{1,3};
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
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_SpreadHigh',num2str(Num)],'-dpng','-r600');
end

viewangle=[65,21;-153,-14;-172,-38;-122,0;81,-13;
    -84,9;-106,36;6,9;-193,-13;9,2];
for Num=1:10
    figure(Num+10)
    SpreadLow=LowData(Num,:);
    SampleNum=SpreadLow{1,1};CellName=SpreadLow{1,2};Frame=SpreadLow{1,3};
    CellIndex=NameDic{find(strcmp(NameDic(:,2),CellName)),1};
    % 形态分割数据
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
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_SpreadLow',num2str(Num)],'-dpng','-r600');
end

