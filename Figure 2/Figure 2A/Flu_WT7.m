% Fluorescence
clear all;clc;close all;

% z layer
layernum=92;
% frame of each embyro
framenumWT=[];
for SampleNum=1:8
    countlist=dir(['..\..\bin\CMap Data\Dataset E\WT_Sample',num2str(SampleNum),'\RawMemb\*.nii.gz']);
    framenumWT=[framenumWT;length(countlist)];
end
TReList=[1.443;1.440;1.443;1.436;1.436;1.403;1.431;1.443];

% WT
figure
SampleNum=7;TRe=TReList(SampleNum);
mkdir(['WT_Sample',num2str(SampleNum)]);
MaxTime=floor(1.44*(framenumWT(SampleNum)-1));
for t=1:framenumWT(SampleNum)
    ImageMemb=zeros(512,712,layernum);
    ImageNuc=zeros(512,712,layernum);
    % Stacking of layers with different focal planes
    for l=1:layernum
        pic=imread(['..\..\bin\CMap Data\Dataset B\WT_Sample',num2str(SampleNum),'\Membrane\WT_Sample',num2str(SampleNum),'-T',num2str(t,'%03d'),'-P',num2str(l,'%02d'),'.tif']);
        ImageMemb(:,:,l)=pic;
        pic=imread(['..\..\bin\CMap Data\Dataset B\WT_Sample',num2str(SampleNum),'\Nucleus\WT_Sample',num2str(SampleNum),'-T',num2str(t,'%03d'),'-P',num2str(l,'%02d'),'.tif']);
        ImageNuc(:,:,l)=pic;
    end
    % corresponding time
    Time=floor(TRe*(t-1));
    Memb=ImageMemb;Nuc=ImageNuc;
    % exclude 0
    FluMembAxes=[];FluNucAxes=[];
    % Traverse the space to find all the serial numbers of the fluorescent spots
    ind=find(Memb>=prctile(reshape(Memb,1,numel(Memb)),95) & Memb>10);[X,Y,Z]=ind2sub(size(Memb),ind);
    FluMembAxes=[X,size(Memb,2)-Y,Z,Memb(ind)];
    ind=find(Nuc>=prctile(reshape(Nuc,1,numel(Nuc)),95) & Nuc>12);[X,Y,Z]=ind2sub(size(Nuc),ind);
    FluNucAxes=[X,size(Nuc,2)-Y,Z,Nuc(ind)];
    % z axis
    Row=find(FluMembAxes(:,3)>=0 & FluMembAxes(:,3)<=layernum);NewMemb=FluMembAxes(Row,:);
    Row=find(FluNucAxes(:,3)>=0 & FluNucAxes(:,3)<=layernum);NewNuc=FluNucAxes(Row,:);
    MaxMemb=max(im2double(NewMemb(:,4)))./4;MaxNuc=max(im2double(NewNuc(:,4)))./2;
    scatter3(NewMemb(:,1),NewMemb(:,2),NewMemb(:,3),im2double(NewMemb(:,4))./MaxMemb,'r','filled','MarkerFaceAlpha',0.02,'MarkerEdgeColor','none');hold on;
    scatter3(NewNuc(:,1),NewNuc(:,2),NewNuc(:,3),im2double(NewNuc(:,4))./MaxNuc,'g','filled','MarkerFaceAlpha',0.02,'MarkerEdgeColor','none');hold on;
    % xlabel('x');ylabel('y');zlabel('z');
    title({['\itt\rm = ',num2str(Time),' min']},'Color','w');
    set(gca,'FontSize',18,'Fontname','Arial');
    axis off;box off;axis equal;
    axis([0 512 0 712 0 92]);
    view([-88,85]);
    set(gcf,'color','k');
    set(gcf,'InvertHardCopy','off');
    set(gca,'position',[0.05 0.05 0.9 0.8]);
    set(gcf,'unit','centimeters','position',[5 8 12 9]);
    % save
    print(gcf,'-dpng',['.\WT_Sample',num2str(SampleNum),'\WT_Sample',num2str(SampleNum),'_',num2str(Time,'%03d'),'.png'],'-r300');
    clf;
end
