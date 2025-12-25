% Neuron: maximum projection sphericity
clear all; clc; close all;

% Name Dictionary
NameDic=readcell(['E:\Project-Gene Expression\CMap Data\Dataset E\name_dictionary.csv']);
NeuronColor=[126,186,210]./255;
PharynxColor=[158,119,221]./255;
MuscleColor=[235,143,182]./255;
IntestineColor=[250,199,86]./255;
SkinColor=[147,202,80]./255;
AxisRotation=[1,-1,-1;1,1,1;1,1,1;-1,-1,1;1,-1,-1;-1,1,-1;1,-1,-1;-1,-1,1];

load(['..\Intestine_ExtremeShapeList.mat']);

% Corey shape factor
HighData=ExtremeFateList{2,1};
LowData=ExtremeFateList{3,1};

viewangle=[20,20;17,10;9,23;0,21;124,-11];
for Num=1:5
    figure(Num)
    % high value
    GSphHigh=HighData(Num,:);
    SampleNum=GSphHigh{1,1};CellName=GSphHigh{1,2};Frame=GSphHigh{1,3};
    CellIndex=NameDic{find(strcmp(NameDic(:,2),CellName)),1};
    load(['..\..\..\bin\Data CMap\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
    Space=AllSpace{1,Frame};
    ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
    CellAxis=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
    [SurfaceAxis,~]=boundary(CellAxis,0.95);
    h=trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',IntestineColor,'FaceAlpha',1);hold on;
    axis off;box off;axis equal;
    view(viewangle(Num,:));
    camlight('headlight','infinite');
    material dull;lighting gouraud;
    title(CellName,'FontSize',15,'FontName','Arial','Color','w');
    % background parameter
    set(gcf,'color','k');
    set(gcf, 'InvertHardCopy','off');
    set(gcf,'unit','centimeters','position',[10,5,8,6]);
    set(gca,'unit','centimeters','position',[2,1,5,3.5]);
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_GSphHigh',num2str(Num)],'-dpng','-r600');
end

viewangle=[168,1;-164,1;-178,16;134,-2;53,40];
for Num=1:5
    figure(Num+5)
    % low value
    GSphLow=LowData(Num,:);
    SampleNum=GSphLow{1,1};CellName=GSphLow{1,2};Frame=GSphLow{1,3};
    CellIndex=NameDic{find(strcmp(NameDic(:,2),CellName)),1};
    load(['..\..\..\bin\CMap Dat\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
    Space=AllSpace{1,Frame};
    ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
    CellAxis=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
    [SurfaceAxis,~]=boundary(CellAxis,0.95);
    h=trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',IntestineColor,'FaceAlpha',1);hold on;
    axis off;box off;axis equal;
    view(viewangle(Num,:));
    camlight('headlight','infinite');
    material dull;lighting gouraud;
    title(CellName,'FontSize',15,'FontName','Arial','Color','w');
    % background parameter
    set(gcf,'color','k');
    set(gcf, 'InvertHardCopy','off');
    set(gcf,'unit','centimeters','position',[10,5,8,6]);
    set(gca,'unit','centimeters','position',[2,1,5,3.5]);
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_GSphLow',num2str(Num)],'-dpng','-r600');
end