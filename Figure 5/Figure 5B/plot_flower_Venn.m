clear all;clc;close all;

load(['.\Fate_GeneMarker.mat']);

XX=GeneMarkerTab(3,[1,2,3,5,4]);

VN=venn(XX{1:5});
% VN=VN.labels(GeneMarkerTab{1,1},GeneMarkerTab{1,2},GeneMarkerTab{1,3},GeneMarkerTab{1,5},GeneMarkerTab{1,4});
VN=VN.draw();

LineList=[47,126,151;
    130,73,143;
    181,59,114;
    229,171,9;
    97,140,44]./255;
FaceList=[185,205,220;
    226,198,232;
    224,190,211;
    245,217,159;
    195,228,174]./255;
for i=1:5
    VN.setPatchN(i,'FaceColor',FaceList(i,:),'EdgeColor',LineList(i,:))
end

set(gca,'FontName','Arial','FontSize',12);
set(gcf,'units','centimeters','Position',[10 5 15 10]);