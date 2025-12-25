% Draw the image of the cell membrane division of the embryo
clear all;clc;close all;

% The rotational view of 18 mechanically-compressed embryos
AxisRotation=[-1,-1,1;
    -1,-1,1;
    1,1,1;
    1,-1,-1;
    -1,-1,1;
    1,-1,-1;
    1,1,1;
    1,-1,-1;
    -1,-1,1;
    -1,1,-1;
    1,1,1;
    1,-1,-1;
    -1,1,-1;
    -1,1,-1;
    -1,1,-1;
    -1,-1,1;
    -1,1,-1];
viewangle=[-90,90;90,-90;-90,180;-90,0;0,180;180,-180];
% name dictionary
NameDic=readcell(['..\..\bin\name_dictionary_CShaper.csv']);
% lineage tree
load('..\..\bin\LineageTree.mat');
% cell color
load(['..\..\bin\ColorTree2.mat']);

for SampleNum=15
    clear AllSpace;
    % data of all frames
    load(['.\Sample',num2str(SampleNum,'%02d'),'_AllSpace.mat']);
    FileNameList={'Left','Right','Dorsal','Ventral','Anterior','Posterior'};
    for i=1:6
        mkdir(['.\CShaper_Sample',num2str(SampleNum),'_',FileNameList{1,i}]);
    end
    FrameNum=length(AllSpace);
    % axis
    meanAxis=[];maxAxis=[];minAxis=[];
    for Frame=1:length(AllSpace)
        Space=AllSpace{1,Frame};
        ExistCell=unique(Space);
        ind=find(Space~=0);[X,Y,Z]=ind2sub(size(Space),ind);
        AllCellAxis=[AxisRotation(SampleNum-3,1)*X,AxisRotation(SampleNum-3,2)*Y,AxisRotation(SampleNum-3,3)*Z];
        % mean/max/min data
        meanAxis=[meanAxis;mean(AllCellAxis(:,1)),mean(AllCellAxis(:,2)),mean(AllCellAxis(:,3))];
        maxAxis=[maxAxis;max(AllCellAxis(:,1)),max(AllCellAxis(:,2)),max(AllCellAxis(:,3))];
        minAxis=[minAxis;min(AllCellAxis(:,1)),min(AllCellAxis(:,2)),min(AllCellAxis(:,3))];
    end
    % centroid
    centroid=mean(meanAxis);
    AllCellAxis=AllCellAxis-centroid;
    CmaxAxis=max(maxAxis)-centroid;CminAxis=min(minAxis)-centroid;
    axisrange=[CminAxis;CmaxAxis]';
    clear AllCellAxis;

    % plot with frame sequencing
    for Frame=1:length(AllSpace)
        Time=round(1.39*(Frame-1));
        Space=AllSpace{1,Frame};
        ExistCell=double(unique(Space));
        ExistCell(find(ExistCell==0))=[];
        CellNumber=length(ExistCell);
        % raw image
        f1=figure(1);
        for Vari=1:length(ExistCell)
            CellIndex=ExistCell(Vari,1);
            % cell identity
            CellName=NameDic{CellIndex+1,2};
            [Row,Column]=find(strcmp(LineageTree,CellName));
            Color=ColorTree{Row,Column};
            ind=find(Space==CellIndex);[X,Y,Z]=ind2sub(size(Space),ind);
            CellAxis=[AxisRotation(SampleNum-3,1)*X,AxisRotation(SampleNum-3,2)*Y,AxisRotation(SampleNum-3,3)*Z];
            CellAxis=CellAxis-centroid;
            [SurfaceAxis,~]=boundary(CellAxis,0.95);
            trisurf(SurfaceAxis,CellAxis(:,1),CellAxis(:,2),CellAxis(:,3),'EdgeColor','none','FaceColor',Color,'FaceAlpha',1);hold on;
            axis off;box off;axis equal;
            axis([axisrange(1,1)-20,axisrange(1,2)+20,axisrange(2,1)-20,axisrange(2,2)+20,axisrange(3,1)-20,axisrange(3,2)+20]);
        end
        % copy image
        for i=2:6
            f2=figure(i);
            copyobj(f1.Children,f2);
        end
        % modify parameter
        for i=1:6
            figure(i)
            view(viewangle(i,:));
            camlight('headlight','infinite');
            material dull;lighting gouraud;            
            title({['\itt\rm = ',num2str(Time),' min'],['Cell Number = ',num2str(CellNumber)]},'FontSize',18,'FontName','Arial','Color','w');
            % background parameter
            set(gcf,'color','k');
            set(gcf, 'InvertHardCopy','off');
            set(gcf,'unit','centimeters','position',[10,5,18,15]);
            set(gca,'unit','centimeters','position',[2,1,14,11]);
            print(gcf,['.\CShaper_Sample',num2str(SampleNum),'_',FileNameList{1,i},'\Sample',num2str(SampleNum),'_Frame_',num2str(Frame,'%03d'),'_Time_',num2str(Time,'%03d')],'-dpng','-r600');
            clf(i);
        end
    end
end
