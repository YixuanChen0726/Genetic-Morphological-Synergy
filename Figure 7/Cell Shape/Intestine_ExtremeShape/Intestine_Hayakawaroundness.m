% Neuron: Hayakawa Roundness
clear all; clc; close all;

% Name Dictionary
NameDic=readcell(['..\..\..\bin\DataSet CMap\name_dictionary.csv']);
NeuronColor=[126,186,210]./255;
PharynxColor=[158,119,221]./255;
MuscleColor=[235,143,182]./255;
IntestineColor=[250,199,86]./255;
SkinColor=[147,202,80]./255;
AxisRotation=[1,-1,-1;1,1,1;1,1,1;-1,-1,1;1,-1,-1;-1,1,-1;1,-1,-1;-1,-1,1];

load(['..\Intestine_ExtremeShapeList.mat']);

% Hayakawa Roundness
HighData=ExtremeFateList{2,5};
LowData=ExtremeFateList{3,5};

viewangle=[-2,28;14,13;9,23;0,21;126,-20];
for Num=1:5
    figure(Num)
    % high value
    HRoundHigh=HighData(Num,:);
    SampleNum=HRoundHigh{1,1};CellName=HRoundHigh{1,2};Frame=HRoundHigh{1,3};
    CellIndex=NameDic{find(strcmp(NameDic(:,2),CellName)),1};
    load(['..\..\..\bin\CMap Data\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
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
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_HRoundHigh',num2str(Num)],'-dpng','-r600');
end

viewangle=[170,2;-181,-9;-67,-14;3,16;-67,-46];
for Num=1:5
    figure(Num+5)
    % low value
    HRoundLow=LowData(Num,:);
    SampleNum=HRoundLow{1,1};CellName=HRoundLow{1,2};Frame=HRoundLow{1,3};
    CellIndex=NameDic{find(strcmp(NameDic(:,2),CellName)),1};
    load(['..\..\..\bin\CMap Data\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
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
    print(gcf,['.\Sample',num2str(SampleNum),'_',CellName,'_Frame_',num2str(Frame,'%03d'),'_HRoundLow',num2str(Num)],'-dpng','-r600');
end
