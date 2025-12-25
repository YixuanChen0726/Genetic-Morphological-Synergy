% Fluorescence
clear all;clc;close all;

SampleNum=15;

% mCherry:Membrane
load(['..\..\bin\CShaper Data\RawMembrane\Sample',num2str(SampleNum,'%02d'),'_AllMemb.mat']);
% gfp:Nucleus
load(['..\..\bin\CShaper Data\RawNucleus\Sample',num2str(SampleNum,'%02d'),'_AllNuc.mat']);

FrameNum=length(AllMemb);MaxTime=floor(1.39*(FrameNum-1));

for Time=0:1:MaxTime
    Frame=round(Time/1.39)+1;
    Memb=AllMemb{1,Frame};Nuc=AllNuc{1,Frame};
    % exclude 0
    FluMembAxes=[];FluNucAxes=[];
    % Traverse the space to find all the serial numbers of the fluorescent spots
    ind=find(Memb>=prctile(reshape(Memb,1,numel(Memb)),93)&Memb>20);[X,Y,Z]=ind2sub(size(Memb),ind);FluMembAxes=[X,Y,Z,Memb(ind)];
    ind=find(Nuc>=prctile(reshape(Nuc,1,numel(Nuc)),93)&Nuc>15);[X,Y,Z]=ind2sub(size(Nuc),ind);FluNucAxes=[X,Y,Z,Nuc(ind)];
    % z axis
    Row=find(FluMembAxes(:,3)>=35 & FluMembAxes(:,3)<=75);NewMemb=FluMembAxes(Row,:);
    Row=find(FluNucAxes(:,3)>=35 & FluNucAxes(:,3)<=75);NewNuc=FluNucAxes(Row,:);
    MaxMemb=max(im2double(NewMemb(:,4)))./10;MaxNuc=max(im2double(NewNuc(:,4)))./5;
    scatter3(NewMemb(:,1),NewMemb(:,2),115-NewMemb(:,3),im2double(NewMemb(:,4))./MaxMemb,'r','filled','MarkerFaceAlpha',0.03,'MarkerEdgeColor','none');hold on;
    scatter3(NewNuc(:,1),NewNuc(:,2),115-NewNuc(:,3),im2double(NewNuc(:,4))./MaxNuc,'g','filled','MarkerFaceAlpha',0.03,'MarkerEdgeColor','none');hold on;
    % xlabel('x');ylabel('y');zlabel('z');
    title({['\itt\rm = ',num2str(Time),' min']},'Color','w');
    set(gca,'FontSize',18,'Fontname','Arial');
    axis off;box off;axis equal;
    axis([0 184 0 256 0 114]);
    view([90,-80]);
    set(gcf,'color','k');
    set(gcf,'InvertHardCopy','off');
    set(gca,'position',[0.05 0.05 0.9 0.8]);
    set(gcf,'unit','centimeters','position',[5 8 12 9]);
    % save data
    print(gcf,'-dpng',['.\Sample',num2str(SampleNum,'%02d'),'\Sample',num2str(SampleNum,'%02d'),'_',num2str(Time,'%03d'),'.png'],'-r600');
    clf;
end
