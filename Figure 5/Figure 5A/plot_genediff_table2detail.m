% The frequency of cell occurrence
clear all;clc;close all;

% cell count
load(['Gene_Diff_Freq.mat']);
GeneDiff=cell2mat(CellDiffTable(2:end,2:end));

% heatmap
FontS=13;
fig=figure('Position',[30,200,1650,600],'Name','CYX');

% Establish the coordinate area
[rows,cols]=size(GeneDiff);

% colorbar region
axBar=subplot(2,3,6,'Position',[1550/1650,50/600,0,500/600]);
axBar.Color='none';
axis off;
uistack(axBar,'bottom');
axBar.Color='none';
axBar.XColor='none';axBar.YColor='none';
CM=colorbar;
CM.Position(3)=0.02;
hold on;
ylabel(CM,'Cell Number','FontSize',FontS,'FontName','Arial');
set(CM,'FontSize',FontS,'FontName','Arial');

% colorbar
c=slanCM('RdPu');
r1=floor(size(c,1)/2);r2=size(c,1)-floor(size(c,1)/2);
c1=interp1([10:1:r1],c(10:r1,:),[10:1:r1]);
c2=interp1([r2:1:size(c,1)],c(r2:size(c,1),:),[r2:0.6:size(c,1)]);
colormap([c1;c(r1:r2-1,:);c2]);
ColorLim=max(max(GeneDiff));
clim([0,ColorLim]);

% gene tree region
axTree2=subplot(2,3,2,'Position',[100/1650,445/600,1300/1650,145/600]);
axTree2.Color='none';
hold on;
% shape tree region
axTree1=subplot(2,3,4,'Position',[10/1650,100/600,85/1650,340/600]);
axTree1.Color='none';
hold on;
% hierarchical clustering
% [observation;variables]
% euclidean distance
ShapeRawDis=pdist(GeneDiff,'euclidean');
GeneRawDis=pdist(GeneDiff.','euclidean');
ShapeDis=squareform(ShapeRawDis);
GeneDis=squareform(GeneRawDis);
ShapeZ=linkage(ShapeDis,'average');
GeneZ=linkage(GeneDis,'average');
% plot clustering tree
rng(0,"twister")
[ShapeFig,~,order1]=dendrogram(ShapeZ,0,'Orientation','left');hold on;
[GeneFig,~,order2]=dendrogram(GeneZ,0,'Orientation','top');hold on;
set(ShapeFig,'Color',[0,0,0],'LineWidth',1.1);
set(GeneFig,'Color',[0,0,0],'LineWidth',1.1);
CellAppear1=GeneDiff(order1,:);
CellAppear2=CellAppear1(:,order2);
shapeName=CellDiffTable(2:end,1);newshapeName=shapeName(order1);
geneName=CellDiffTable(1,2:end);newgeneName=geneName(order2);

% modification
tempFig1=ShapeFig(1).Parent.Parent;
axTree1=copyAxes(tempFig1,1,axTree1);
axTree1.XColor='none';axTree1.YColor='none';
axTree1.YDir='reverse';
axTree1.XTick=[];axTree1.YTick=[];
% Ylimitation
b=find(strcmp(newshapeName,'Neuron / Pharynx'));
axTree1.YLim=[b,rows]+[-.5,.5];
delete(tempFig1);
tempFig2=GeneFig(1).Parent.Parent;
axTree2=copyAxes(tempFig2,1,axTree2);
axTree2.XColor='none';axTree2.YColor='none';
axTree2.XTick=[];axTree2.YTick=[];
a=find(strcmp(newgeneName,'aly-2'));
axTree2.XLim=[a,cols]+[-.5,.5];
delete(tempFig2);

% heatmap region
axMain=subplot(2,3,5,'Position',[100/1650,100/600,1300/1650,340/600]);
axMain.Color='none';
axMain.YDir='reverse';
% heatmap
for i=b:1:size(CellAppear2,1)
    for j=a:1:size(CellAppear2,2)
        patch(axMain,[j-0.5,j-0.5,j+0.5,j+0.5],[i-0.5,i+0.5,i+0.5,i-0.5],CellAppear2(i,j),'EdgeColor','none');hold on;
    end
end

xticks([a:1:size(CellAppear2,2)]);xtickangle(90);
axMain.XTickLabel=arrayfun(@(x) ['\it ',x{:}],newgeneName(1,a:end),'UniformOutput',false);
yticks([b:1:size(CellAppear2,1)]);ytickangle(0);
axMain.YAxisLocation='right';
set(axMain,'YTickLabel',newshapeName(b:end,1),'FontAngle','normal');
axis([a-0.5,size(GeneDiff,2)+0.5,b-0.5,size(GeneDiff,1)+0.5]);
set(axMain,'FontName','Arial','FontSize',FontS);

function axbag=copyAxes(fig,k,newAx)
classList(length(fig.Children))=true;
for n=1:length(fig.Children)
    classList(n)=isa(fig.Children(n),'matlab.graphics.axis.Axes');
end
isaaxes=find(classList);
oriAx=fig.Children(isaaxes(end-k+1));
if isaaxes(end-k+1)-1<1 || isa(fig.Children(isaaxes(end-k+1)-1),'matlab.graphics.axis.Axes')
    oriLgd=[];
else
    oriLgd=fig.Children(isaaxes(end-k+1)-1);
end
axbag=copyobj([oriAx,oriLgd],newAx.Parent);
axbag(1).Position=newAx.Position;
delete(newAx)
end

function colorList=slanCM(type,num)
if nargin<2
    num=256;
end
if nargin<1
    type='';
end

slanCM_Data=load('..\..\Figure 2\slanCM_Data.mat');
CList_Data=[slanCM_Data.slandarerCM(:).Colors];

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
