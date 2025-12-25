% CMap v.s. CShaper : each point for one cell
clear all;clc;close all;

load('ABplpapp_All_ContactRelation.mat');
AllContact=AllContact([1:26],:);
% delete zero
Allzero=[];
for i=2:size(AllContact,2)
    data=cell2mat(AllContact(2:end,i));
    if all(data(:) == 0)
        Allzero=[Allzero,i];
    end
end
AllContact(:,Allzero)=[];

CMapContact=AllContact([1:9],:);
CShaperContact=AllContact([1,10:26],:);

% 在下方计算均值
CMapContact{size(CMapContact,1)+1,1}='mean contact area';
for i=2:size(CMapContact,2)
    data=cell2mat(CMapContact(2:end,i));
    CMapContact{end,i}=mean(data);
end

CShaperContact{size(CShaperContact,1)+1,1}='mean contact area';
for i=2:size(CShaperContact,2)
    data=cell2mat(CShaperContact(2:end,i));
    CShaperContact{end,i}=mean(data);
end

% data
CMapdata=cell2mat(CMapContact(size(CMapContact,1),2:end));
CShaperdata=cell2mat(CShaperContact(size(CShaperContact,1),2:end));
% MSapp,MSappp,MSappa,MSapap
MSRow=[find(contains(CMapContact(1,2:end),'MS'))];
% plot
Grey=[150,150,150]./255;
Red=[221,75,74]./255;
x=-10:0.1:40;y=x;
plot(x,y,'--','Color','k','LineWidth',1.5);hold on;
for i=1:length(CMapdata)
    if ~isempty(intersect(i,MSRow))
        scatter(CMapdata(i),CShaperdata(i),35,Red,'filled','MarkerFaceAlpha',0.7);hold on;
    else
        scatter(CMapdata(i),CShaperdata(i),35,Grey,'filled','MarkerFaceAlpha',0.6);hold on;
    end
end

% goodness of fit
SSE=sum((CShaperdata-CMapdata).^2);
SST=sum((CShaperdata-mean(CShaperdata)).^2);
RSquare=1-SSE/SST;
text(0.1,27,['\itG\rm = ',num2str(RSquare,'%.3f')],'FontName','Arial','FontSize',13);

xticks(0:10:30);yticks(0:10:30);
xlabel({'Cell-cell contact area under','natural condition (\mum^2)'});
ylabel({'Cell-cell contact area with ','mechanical compression (\mum^2)'});
axis([-2,30,-2,30]);
axis square;
set(gcf,'units','centimeters','position',[10,5,12,10]);
set(gca,'units','centimeters','position',[1.5,2.5,10,7]);
set(gca,'FontName','Arial','FontSize',13);
