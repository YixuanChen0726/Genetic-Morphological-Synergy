% the membrane segmentation of embryo with a white background
clear all;clc;close all;

% rotation relation of 8 embryos
AxisRotation=[1,-1,-1;1,1,1;1,1,1;-1,-1,1;1,-1,-1;-1,1,-1;1,-1,-1;-1,-1,1];
viewangle=[-90,90;90,-90;-90,180;-90,0;0,180;180,-180];
% name dictionary
NameDic=readcell(['..\..\bin\DataSet CMap\name_dictionary.csv']);
% lineage tree
load('..\..\bin\LineageTree.mat');
% fate distribution of all cells
load('..\..\bin\FateClassification.mat');

NeuronColor=[106,195,202]./255;
PharynxColor=[158,119,221]./255;
MuscleColor=[235,143,182]./255;
IntestineColor=[250,199,86]./255;
SkinColor=[147,202,80]./255;
GermColor=[103,159,227]./255;
DeathColor=[215,111,111]./255;
OtherColor=[200,200,200]./255;

FateColor=[NeuronColor;PharynxColor;MuscleColor;IntestineColor;SkinColor;GermColor;DeathColor;OtherColor];
FateNameList={'Neuron','Pharynx','Muscle','Intestine','Skin','Germ Cell','Death','Unspecified/Other'};

Mothername='EMS';
dau1name=  'MS';
dau2name=  'E';

% position of three cells
[MRow,MColumn]=find(strcmp(LineageTree,Mothername));
[D1Row,D1Column]=find(strcmp(LineageTree,dau1name));
[D2Row,D2Column]=find(strcmp(LineageTree,dau2name));

% change Daughter1 sublineage cell list to cell index
Dau1List={};
for DauLayer=D1Column:1:size(LineageTree,2)
    n=2^(DauLayer-D1Column);
    for DauIndex=(n*D1Row-(n-1)):1:min([(n*D1Row),size(LineageTree,1)])
        DauName=LineageTree{DauIndex,DauLayer};
        if all(DauName)
            Dau1List=[Dau1List;{DauName}];

        end
    end
end
% MS:pharynx fate
MSPharynx={};
Row1=find(strcmp(FateClas(:,2),'Pharynx'));
PharynxList=FateClas(Row1,:);
for i=1:length(Dau1List)
    NameTemp=Dau1List{i,1};
    k=find(strcmp(PharynxList(:,1),NameTemp));
    if ~isempty(k)
        MSPharynx=[MSPharynx;PharynxList(k,1)];
    end
end
Dau1PharynxIndex=[];
for i=1:length(MSPharynx)
    CellName=MSPharynx{i,1};
    Row=find(strcmp(NameDic(:,2),CellName));
    Dau1PharynxIndex=[Dau1PharynxIndex;NameDic{Row,1}];
end
% MS:muscle fate
MSMuscle={};
Row1=find(strcmp(FateClas(:,2),'Muscle'));
MuscleList=FateClas(Row1,:);
for i=1:length(Dau1List)
    NameTemp=Dau1List{i,1};
    k=find(strcmp(MuscleList(:,1),NameTemp));
    if ~isempty(k)
        MSMuscle=[MSMuscle;MuscleList(k,1)];
    end
end
Dau1MuscleIndex=[];
for i=1:length(MSMuscle)
    CellName=MSMuscle{i,1};
    Row=find(strcmp(NameDic(:,2),CellName));
    Dau1MuscleIndex=[Dau1MuscleIndex;NameDic{Row,1}];
end

% change Daughter2 sublineage cell list to cell index
Dau2List={};
for DauLayer=D2Column:1:size(LineageTree,2)
    n=2^(DauLayer-D2Column);
    for DauIndex=(n*D2Row-(n-1)):1:min([(n*D2Row),size(LineageTree,1)])
        DauName=LineageTree{DauIndex,DauLayer};
        if all(DauName)
            Dau2List=[Dau2List;{DauName}];
        end
    end
end
% E:intestine fate
EIntestine={};
Row1=find(strcmp(FateClas(:,2),'Intestine'));
NeuronList=FateClas(Row1,:);
for i=1:length(Dau2List)
    NameTemp=Dau2List{i,1};
    k=find(strcmp(NeuronList(:,1),NameTemp));
    if ~isempty(k)
        EIntestine=[EIntestine;NeuronList(k,1)];
    end
end
Dau2Index=[];
for i=1:length(EIntestine)
    CellName=EIntestine{i,1};
    Row=find(strcmp(NameDic(:,2),CellName));
    Dau2Index=[Dau2Index;NameDic{Row,1}];
end

figure
for SampleNum=1:8
    load(['..\..\bin\CMap Data\Dataset E\WT_Sample',num2str(SampleNum),'_AllSpace.mat']);
    mkdir(['.\Sample',num2str(SampleNum)]);
    FrameNum=length(AllSpace);MaxTime=floor(1.44*(FrameNum-1));
    % axis range
    meanAxis=[];maxAxis=[];minAxis=[];
    for Frame=1:length(AllSpace)
        Space=AllSpace{1,Frame};
        ExistCell=unique(Space);
        ExistCell(find(ExistCell==0))=[];
        AllCellAxis=[];
        for Vari=1:length(ExistCell)
            CellIndex=ExistCell(Vari);
            ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
            AllCellAxis=[AllCellAxis;AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
        end
        % mean,max,min
        meanAxis=[meanAxis;mean(AllCellAxis(:,1)),mean(AllCellAxis(:,2)),mean(AllCellAxis(:,3))];
        maxAxis=[maxAxis;max(AllCellAxis(:,1)),max(AllCellAxis(:,2)),max(AllCellAxis(:,3))];
        minAxis=[minAxis;min(AllCellAxis(:,1)),min(AllCellAxis(:,2)),min(AllCellAxis(:,3))];
    end
    % centroid
    centroid=mean(meanAxis,1);
    CmaxAxis=max(maxAxis,[],1)-centroid;CminAxis=min(minAxis,[],1)-centroid;
    axisrange=[CminAxis;CmaxAxis]';

    % plot
    for Frame=1:length(AllSpace)
        Space=AllSpace{1,Frame};
        ExistCell=double(unique(Space));
        ExistCell(find(ExistCell==0))=[];
        CellNumber=length(ExistCell)-1;
        DauNum=[Dau1PharynxIndex;Dau1MuscleIndex;Dau2Index];
        if ~isempty(intersect(ExistCell,DauNum))
            clf;
            set(gcf,'unit','centimeters','position',[2 8 40 12]);
            t=tiledlayout(1,7);
            t.TileSpacing = 'compact';%由宽至窄 'loose'|'compact'|'tight'|'none'
            t.Padding = 'compact';%由大面积填充到紧致填充 'loose'|'compact'|'tight'
            titlestr={'A-P/D-V(L)','A-P/D-V(R)','A-P/L-R(D)','A-P/L-R(V)','D-V/L-R(A)','D-V/L-R(P)'};
            % the first figure
            ax=nexttile(1);
            for Vari=1:length(ExistCell)
                CellIndex=ExistCell(Vari,1);
                if ~isempty(intersect(CellIndex,DauNum))
                    alpha=1;
                else
                    alpha=0.05;
                end
                CellName=NameDic{CellIndex+1,2};Color=OtherColor; 
                if ~isempty(find(strcmp(FateClas(:,1),CellName)))
                    CellFate=FateClas{find(strcmp(FateClas(:,1),CellName)),2};
                    Row=find(strcmp(FateNameList,CellFate));
                    if ~isempty(Row)
                        Color=FateColor(Row,:);
                    else
                        Color=OtherColor;
                    end
                end
                ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
                CellAxis=[AxisRotation(SampleNum,1)*X,AxisRotation(SampleNum,2)*Y,AxisRotation(SampleNum,3)*Z];
                CellAxis=CellAxis-centroid;
                [SurfaceAxis,~]=boundary(CellAxis,0.95);
                h=trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',Color,'FaceAlpha',alpha);hold on;
            end
            axis off;box off;axis equal;
            axis([axisrange(1,1)-20,axisrange(1,2)+20,axisrange(2,1)-20,axisrange(2,2)+20,axisrange(3,1)-20,axisrange(3,2)+20]);
            for i=2:6
                axc=copyobj(ax,t);
                ax(1).Layout.Tile=i;
            end
            % modify image parameter
            for i=1:6
                nexttile(i);
                material dull;lighting gouraud;
                view(viewangle(i,:));
                camlight('headlight','infinite');
                if strcmp(titlestr{1,i},'A-P/D-V(L)')
                    text('string',Mothername,'Units','normalized','position',[-0.1,0.2],'FontSize',18,'FontName','Arial','Color','k','Rotation',90);
                end
            end
            % legend
            nexttile(7);
            hx=fliplr(1:1:12);
            for k=1:6
                patch([0 2 2 0],[hx(2*k-1),hx(2*k-1),hx(2*k),hx(2*k)],FateColor(k,:),'EdgeColor','none');hold on;
                text(2.5,(hx(2*k-1)+hx(2*k))./2,[FateNameList{1,k}],'FontSize',15,'FontName','Arial','Color','k');
            end
            axis([-0.2,13,-1,13]);axis off;
            title(t,{['Cell Number = ',num2str(CellNumber)],[],['  ',titlestr{1,1},'              ',titlestr{1,2},'              ',titlestr{1,3},'              ',titlestr{1,4},'              ',titlestr{1,5},'              ',titlestr{1,6},'                              ']},'FontSize',18,'FontName','Arial','Color','k');
            % background parameter
            set(gcf,'color','w');
            set(gcf, 'InvertHardCopy','off');
            set(t,'unit','centimeters','position',[2,1,37,7]);
            print(gcf,['.\Sample',num2str(SampleNum),'\Sample',num2str(SampleNum),'_',num2str(Frame,'%03d')],'-dpng','-r600');
            clf;
        end
    end
end
