% Ring-shaped volcanic map
% comparison of pharynx, muscle and intestine fate of EMS
clear all;clc; close all;

% EMS
% asymmetry of all cells of sublineages
load(['..\..\bin\Coasymmetry\all gene\EMS_AllGeneSig.mat']);
EMSAll=AllGeneSigni;
% asymmetry of partial cells of specific fates
% Pharynx v.s. Intestine
load(['..\..\bin\Coasymmetry\partial fate gene\EMS_Pharynx_Intestine.mat']);
EMSPartial1=PartialGeneSigni;
% Muscle v.s. Intestine
load(['..\..\bin\Coasymmetry\partial fate gene\EMS_Muscle_Intestine.mat']);
EMSPartial2=PartialGeneSigni;

% Classification
GreyData1={};LinG1=0;
GreyData2={};LinG2=0;
% Pharynx v.s. Intestine
AnEmphData1={};LinAE1=0;PoEmphData1={};LinPE1=0;
AnSinData1={};LinAS1=0;PoSinData1={};LinPS1=0;
TextData1={};LinT1=0;
% Muscle v.s. Intestine
AnEmphData2={};LinAE2=0;PoEmphData2={};LinPE2=0;
AnSinData2={};LinAS2=0;PoSinData2={};LinPS2=0;

for i=1:size(AllGeneSigni,2)
    % shape name
    GeneName=EMSPartial1{1,i};
    % mean data
    MeanGeneAll=EMSAll{4,i};
    MeanGenePar1=EMSPartial1{4,i};
    MeanGenePar2=EMSPartial2{4,i};
    % p value and asymmetry of all and partial cells
    PsigPar1=EMSPartial1{5,i};
    AsymPar1=EMSPartial1{6,i};
    PsigPar2=EMSPartial2{5,i};
    AsymPar2=EMSPartial2{6,i};
    PsigAll=EMSAll{5,i};
    AsymAll=EMSAll{6,i};
    % Pharynx v.s. Intestine
    % p value (all & partial) < 0.1
    if (PsigPar1<0.1) && (PsigAll<0.1)
        % p value decrease and asymmetry increase
        if (PsigPar1<PsigAll) && (abs(AsymPar1)>=abs(AsymAll))
            % Anterior:Pharynx highly expressed
            if AsymPar1>0
                LinAE1=LinAE1+1;
                AnEmphData1(:,LinAE1)={GeneName;AsymPar1;PsigPar1;MeanGenePar1(1,1);MeanGenePar1(2,1)};
            end
            % Posterior:Intestine highly expressed
            if AsymPar1<0
                LinPE1=LinPE1+1;
                PoEmphData1(:,LinPE1)={GeneName;AsymPar1;PsigPar1;MeanGenePar1(1,1);MeanGenePar1(2,1)};
            end
        % p value decrease or asymmetry increase
        elseif (PsigPar1<PsigAll) || (abs(AsymPar1)>=abs(AsymAll))
            % Anterior:Pharynx highly expressed
            if AsymPar1>0
                LinAS1=LinAS1+1;
                AnSinData1(:,LinAS1)={GeneName;AsymPar1;PsigPar1;MeanGenePar1(1,1);MeanGenePar1(2,1)};
            end
            % Posterior:Intestine highly expressed
            if AsymPar1<0
                LinPS1=LinPS1+1;
                PoSinData1(:,LinPS1)={GeneName;AsymPar1;PsigPar1;MeanGenePar1(1,1);MeanGenePar1(2,1)};
            end
        % p value increase and asymmetry decrease
        else
            LinG1=LinG1+1;
            GreyData1(:,LinG1)={GeneName;AsymPar1;PsigPar1;MeanGenePar1(1,1);MeanGenePar1(2,1)};
        end
    else
        LinG1=LinG1+1;
        GreyData1(:,LinG1)={GeneName;AsymPar1;PsigPar1;MeanGenePar1(1,1);MeanGenePar1(2,1)};
    end
    % TextList
    if (PsigPar1<1e-40) && (PsigAll<1e-3) && (abs(AsymPar1)>=0.3)
        LinT1=LinT1+1;
        TextData1(:,LinT1)={GeneName;AsymPar1;PsigPar1;MeanGenePar1(1,1);MeanGenePar1(2,1)};
    end

    % Muscle v.s. Intestine
    % p value (all & partial) < 0.1
    if (PsigPar2<0.1) && (PsigAll<0.1)
        % p value decrease and asymmetry increase
        if (PsigPar2<PsigAll) && (abs(AsymPar2)>=abs(AsymAll))
            % Anterior:Muscle highly expressed
            if AsymPar2>0
                LinAE2=LinAE2+1;
                AnEmphData2(:,LinAE2)={GeneName;AsymPar2;PsigPar2;MeanGenePar2(1,1);MeanGenePar2(2,1)};
            end
            % Posterior:Intestine highly expressed
            if AsymPar2<0
                LinPE2=LinPE2+1;
                PoEmphData2(:,LinPE2)={GeneName;AsymPar2;PsigPar2;MeanGenePar2(1,1);MeanGenePar2(2,1)};
            end
         % p value decrease or asymmetry increase
        elseif (PsigPar2<PsigAll) || (abs(AsymPar2)>=abs(AsymAll))
            % Anterior:Muscle highly expressed
            if AsymPar2>0
                LinAS2=LinAS2+1;
                AnSinData2(:,LinAS2)={GeneName;AsymPar2;PsigPar2;MeanGenePar2(1,1);MeanGenePar2(2,1)};
            end
            % Posterior:Intestine highly expressed
            if AsymPar2<0
                LinPS2=LinPS2+1;
                PoSinData2(:,LinPS2)={GeneName;AsymPar2;PsigPar2;MeanGenePar2(1,1);MeanGenePar2(2,1)};
            end
        % p value increase and asymmetry decrease
        else
            LinG2=LinG2+1;
            GreyData2(:,LinG2)={GeneName;AsymPar2;PsigPar2;MeanGenePar2(1,1);MeanGenePar2(2,1)};
        end
    else
        LinG2=LinG2+1;
        GreyData2(:,LinG2)={GeneName;AsymPar2;PsigPar2;MeanGenePar2(1,1);MeanGenePar2(2,1)};
    end
end

% plot
LinW=1.5;MarS=33;Alp=0.6;dx=45;
% double weaken
GreyBG=[180,180,180]./255;
% double enhance
IntestineColor=[229,171,9]./255;
PharynxColor=[158,119,221]./255;
MuscleColor=[235,143,182]./255;
% single enhance
SInteColor=[250,199,86]./255;
SPhaColor=[226,198,232]./255;
SMusColor=[244,200,221]./255;

% x:asymmetry; y:p value
% double weaken:Pharynx
AsymValue=cell2mat(GreyData1(2,:));
PValue=-log10(cell2mat(GreyData1(3,:)))+dx;
scatter(AsymValue,PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',GreyBG,'MarkerEdgeColor',GreyBG);hold on;
% Muscle
AsymValue=cell2mat(GreyData2(2,:));
% p value
PValue=-log10(cell2mat(GreyData2(3,:)))+dx;
scatter(AsymValue,-PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',GreyBG,'MarkerEdgeColor',GreyBG);hold on;

% single strength:Pharynx
AsymValue=cell2mat(AnSinData1(2,:));
PValue=-log10(cell2mat(AnSinData1(3,:)))+dx;
scatter(AsymValue,PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',PharynxColor,'MarkerEdgeColor',PharynxColor);hold on;
AsymValue=cell2mat(PoSinData1(2,:));
PValue=-log10(cell2mat(PoSinData1(3,:)))+dx;
scatter(AsymValue,PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',IntestineColor,'MarkerEdgeColor',IntestineColor);hold on;

% single strength:Muscle
AsymValue=cell2mat(AnSinData2(2,:));
PValue=-log10(cell2mat(AnSinData2(3,:)))+dx;
scatter(AsymValue,-PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',MuscleColor,'MarkerEdgeColor',MuscleColor);hold on;
AsymValue=cell2mat(PoSinData2(2,:));
PValue=-log10(cell2mat(PoSinData2(3,:)))+dx;
scatter(AsymValue,-PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',IntestineColor,'MarkerEdgeColor',IntestineColor);hold on;

% double strength:Pharynx
AsymValue=cell2mat(AnEmphData1(2,:));
PValue=-log10(cell2mat(AnEmphData1(3,:)))+dx;
scatter(AsymValue,PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',PharynxColor,'MarkerEdgeColor',PharynxColor);hold on;
AsymValue=cell2mat(PoEmphData1(2,:));
PValue=-log10(cell2mat(PoEmphData1(3,:)))+dx;
scatter(AsymValue,PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',IntestineColor,'MarkerEdgeColor',IntestineColor);hold on;

% double strength:Muscle
AsymValue=cell2mat(AnEmphData2(2,:));
PValue=-log10(cell2mat(AnEmphData2(3,:)))+dx;
scatter(AsymValue,-PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',MuscleColor,'MarkerEdgeColor',MuscleColor);hold on;
AsymValue=cell2mat(PoEmphData2(2,:));
PValue=-log10(cell2mat(PoEmphData2(3,:)))+dx;
scatter(AsymValue,-PValue,MarS,'filled','MarkerFaceAlpha',Alp,'MarkerFaceColor',IntestineColor,'MarkerEdgeColor',IntestineColor);hold on;

% subtitle
text(-0.65,290+dx,{'Intestine highly expressed','Pharynx low expressed'},'FontName','Arial','FontSize',12,'HorizontalAlignment','center');
text(0.65,290+dx,{'Pharynx highly expressed','Intestine low expressed'},'FontName','Arial','FontSize',12,'HorizontalAlignment','center');
text(-0.65,-290-dx,{'Intestine highly expressed','Muscle low expressed'},'FontName','Arial','FontSize',12,'HorizontalAlignment','center');
text(0.65,-290-dx,{'Muscle highly expressed','Intestine low expressed'},'FontName','Arial','FontSize',12,'HorizontalAlignment','center');

ptheshold=-log10(0.1)+dx;Atheshold=0.1;
plot([-5 5],[ptheshold,ptheshold],'--','Color','k','LineWidth',LinW);hold on;
plot([-5 5],[-ptheshold,-ptheshold],'--','Color','k','LineWidth',LinW);hold on;
plot([-Atheshold -Atheshold],[-400-dx,400+dx],'--','Color','k','LineWidth',LinW);hold on;
plot([Atheshold Atheshold],[-400-dx,400+dx],'--','Color','k','LineWidth',LinW);hold on;
xlabel('(\itG\rm_{MS}-\itG\rm_E)/(\itG\rm_{MS}+\itG\rm_E)');
yticks([[-250:50:0]-dx,[0:50:250]+dx]);
yticklabels({'250','200','150','100','50','0','0','50','100','150','200','250'});
ylabel('-log_{10}(\itp\rm-value)        ');
axis([-1.2,1.2,-360,360]);
box on;

set(gca,'XAxisLocation','origin');
set(gca,'FontName','Arial','FontSize',12);
set(gca,'units','centimeters','position',[2.2 1 12.5 22.5]);
set(gcf,'units','centimeters','position',[10 2 15 24]);

