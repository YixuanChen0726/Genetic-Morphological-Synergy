% EMS shape comparison
clear all; clc; close all;

% all cell asymmetry
load(['..\..\bin\Coasymmetry\all shape\EMS_AllShapeSig.mat']);
EMSAll=AllShapeSigni;
% Pharynx asymmetry
load(['..\..\bin\Coasymmetry\partial fate shape\EMS_Pharynx_Intestine.mat']);
EMSPartial1=PartialShapeSigni;
% Muscle asymmetry
load(['..\..\bin\Coasymmetry\partial fate shape\EMS_Muscle_Intestine.mat']);
EMSPartial2=PartialShapeSigni;

PlotData={};LineIndex=0;
for i=1:size(EMSAll,2)
    % shape name
    ShapeName=EMSPartial1{1,i};
    % p value and asymmetry of all cells
    PsigPar1=EMSPartial1{5,i};
    AsymPar1=EMSPartial1{6,i};
    PsigPar2=EMSPartial2{5,i};
    AsymPar2=EMSPartial2{6,i};
    PsigAll=EMSAll{5,i};
    AsymAll=EMSAll{6,i};
    % calculate partial mean/std
    PharynxData=EMSPartial1{2,i};
    MuscleData=EMSPartial2{2,i};
    IntestineData=EMSPartial1{3,i};
    % scale
    AllData=mapminmax([cell2mat(PharynxData(:,2))',cell2mat(MuscleData(:,2))',cell2mat(IntestineData(:,2))'],0,1);
    ScPharynxData=AllData(1,1:size(PharynxData,1));
    ScMuscleData=AllData(1,size(PharynxData,1)+1:size(PharynxData,1)+size(MuscleData,1));
    ScIntestineData=AllData(1,size(PharynxData,1)+size(MuscleData,1)+1:end);
    % histogram:pharynx,intestine, muscle
    meandata=[mean(ScPharynxData);mean(ScIntestineData);mean(ScMuscleData)];
    stddata=[std(ScPharynxData);std(ScIntestineData);std(ScMuscleData)];
    % p value (all & partial) < 0.1
    if ((PsigPar1<0.1) && (PsigPar2<0.1)) && (PsigAll<0.1)
        % p value decrease && asymmetry increase
        if ((PsigPar1<PsigAll) || (abs(AsymPar1)>=abs(AsymAll))) || ((PsigPar2<PsigAll) || (abs(AsymPar2)>=abs(AsymAll)))
            LineIndex=LineIndex+1;
            PlotData{1,LineIndex}=ShapeName;
            % partial p value
            PlotData{2,LineIndex}=[PsigPar1;PsigPar2];
            % Q1
            PlotData{3,LineIndex}=meandata;
            % Q2
            PlotData{4,LineIndex}=stddata;
        end
    end
end

% circle boxplot
LinW=1.1;FSize=15;
PharynxColor=[158,119,221]./255;
MuscleColor=[235,143,182]./255;
IntestineColor=[250,199,86]./255;

r=0.3;
for i=1:size(PlotData,2)
    % angle
    dangle=2*pi/size(PlotData,2);
    % Pharynx:mean
    theta1=linspace((i-1)*dangle,(i-1)*dangle+dangle/4,100);
    m1=PlotData{3,i}(1,1);s1=PlotData{4,i}(1,1);
    fill([r*cos(theta1),fliplr((r+m1)*cos(theta1))],[r*sin(theta1),fliplr((r+m1)*sin(theta1))],PharynxColor,'EdgeColor','k','LineWidth',LinW,'FaceAlpha',0.7);hold on;
    % Pharynx:std
    beta1=(i-1)*dangle+dangle/8;dt=0.03;
    plot([(r+m1-s1/2)*cos(beta1),(r+m1+s1/2)*cos(beta1)],[(r+m1-s1/2)*sin(beta1),(r+m1+s1/2)*sin(beta1)],'Color','k','LineWidth',LinW);hold on;
    plot([(r+m1-s1/2)*cos(beta1-dt),(r+m1-s1/2)*cos(beta1+dt)],[(r+m1-s1/2)*sin(beta1-dt),(r+m1-s1/2)*sin(beta1+dt)],'Color','k','LineWidth',LinW);hold on;
    plot([(r+m1+s1/2)*cos(beta1-dt),(r+m1+s1/2)*cos(beta1+dt)],[(r+m1+s1/2)*sin(beta1-dt),(r+m1+s1/2)*sin(beta1+dt)],'Color','k','LineWidth',LinW);hold on;
    % Intestine:mean
    theta2=linspace((i-1)*dangle+dangle/4,(i-1)*dangle+dangle/2,100);
    m2=PlotData{3,i}(2,1);s2=PlotData{4,i}(2,1);
    fill([r*cos(theta2),fliplr((r+m2)*cos(theta2))],[r*sin(theta2),fliplr((r+m2)*sin(theta2))],IntestineColor,'EdgeColor','k','LineWidth',LinW,'FaceAlpha',0.7);hold on;
    % Intestine:std
    beta2=(i-1)*dangle+3*dangle/8;
    plot([(r+m2-s2/2)*cos(beta2),(r+m2+s2/2)*cos(beta2)],[(r+m2-s2/2)*sin(beta2),(r+m2+s2/2)*sin(beta2)],'Color','k','LineWidth',LinW);hold on;
    plot([(r+m2-s2/2)*cos(beta2-dt),(r+m2-s2/2)*cos(beta2+dt)],[(r+m2-s2/2)*sin(beta2-dt),(r+m2-s2/2)*sin(beta2+dt)],'Color','k','LineWidth',LinW);hold on;
    plot([(r+m2+s2/2)*cos(beta2-dt),(r+m2+s2/2)*cos(beta2+dt)],[(r+m2+s2/2)*sin(beta2-dt),(r+m2+s2/2)*sin(beta2+dt)],'Color','k','LineWidth',LinW);hold on;
    % Muscle:mean
    theta3=linspace((i-1)*dangle+dangle/2,(i-1)*dangle+3*dangle/4,100);
    m3=PlotData{3,i}(3,1);s3=PlotData{4,i}(3,1);
    fill([r*cos(theta3),fliplr((r+m3)*cos(theta3))],[r*sin(theta3),fliplr((r+m3)*sin(theta3))],MuscleColor,'EdgeColor','k','LineWidth',LinW,'FaceAlpha',0.7);hold on;
    % Muscle:std
    beta3=(i-1)*dangle+5*dangle/8;
    plot([(r+m3-s3/2)*cos(beta3),(r+m3+s3/2)*cos(beta3)],[(r+m3-s3/2)*sin(beta3),(r+m3+s3/2)*sin(beta3)],'Color','k','LineWidth',LinW);hold on;
    plot([(r+m3-s3/2)*cos(beta3-dt),(r+m3-s3/2)*cos(beta3+dt)],[(r+m3-s3/2)*sin(beta3-dt),(r+m3-s3/2)*sin(beta3+dt)],'Color','k','LineWidth',LinW);hold on;
    plot([(r+m3+s3/2)*cos(beta3-dt),(r+m3+s3/2)*cos(beta3+dt)],[(r+m3+s3/2)*sin(beta3-dt),(r+m3+s3/2)*sin(beta3+dt)],'Color','k','LineWidth',LinW);hold on;
    % p value comparison
    gamma1=linspace(beta1-dt,beta2-dt,20);
    gamma2=linspace(beta2+dt,beta3+dt,20);
    % Pharynx & Intestine
    plot((r+0.8)*cos(gamma1),(r+0.8)*sin(gamma1),'Color','k','LineWidth',LinW);hold on;
    plot([(r+0.75)*cos(beta1-dt),(r+0.8)*cos(beta1-dt)],[(r+0.75)*sin(beta1-dt),(r+0.8)*sin(beta1-dt)],'Color','k','LineWidth',LinW);hold on;
    plot([(r+0.75)*cos(beta2-dt),(r+0.8)*cos(beta2-dt)],[(r+0.75)*sin(beta2-dt),(r+0.8)*sin(beta2-dt)],'Color','k','LineWidth',LinW);hold on;
    % Muscle & Intestine
    plot((r+0.8)*cos(gamma2),(r+0.8)*sin(gamma2),'Color','k','LineWidth',LinW);hold on;
    plot([(r+0.75)*cos(beta3+dt),(r+0.8)*cos(beta3+dt)],[(r+0.75)*sin(beta3+dt),(r+0.8)*sin(beta3+dt)],'Color','k','LineWidth',LinW);hold on;
    plot([(r+0.75)*cos(beta2+dt),(r+0.8)*cos(beta2+dt)],[(r+0.75)*sin(beta2+dt),(r+0.8)*sin(beta2+dt)],'Color','k','LineWidth',LinW);hold on;
    % p value
    ppharynx=PlotData{2,i}(1,1);rho1=(i-1)*dangle+dangle/4-dt;
    if ppharynx<=1e-4
        text((r+0.86)*cos(rho1),(r+0.86)*sin(rho1),'****','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho1*(180/pi));
    elseif ppharynx>=1e-4 && ppharynx<1e-3
        text((r+0.86)*cos(rho1),(r+0.86)*sin(rho1),'***','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho1*(180/pi));
    elseif ppharynx>=1e-3 && ppharynx<1e-2
        text((r+0.86)*cos(rho1),(r+0.86)*sin(rho1),'**','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho1*(180/pi));
    elseif ppharynx>=1e-2 && ppharynx<1e-1
        text((r+0.86)*cos(rho1),(r+0.86)*sin(rho1),'*','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho1*(180/pi));
    else
        text((r+0.86)*cos(rho1),(r+0.86)*sin(rho1),'n.s.','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho1*(180/pi));
    end
    pmuscle=PlotData{2,i}(2,1);rho2=(i-1)*dangle+dangle/2+dt;
    if pmuscle<=1e-4
        text((r+0.86)*cos(rho2),(r+0.86)*sin(rho2),'****','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho2*(180/pi));
    elseif pmuscle>=1e-4 && pmuscle<1e-3
        text((r+0.86)*cos(rho2),(r+0.86)*sin(rho2),'***','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho2*(180/pi));
    elseif pmuscle>=1e-3 && pmuscle<1e-2
        text((r+0.86)*cos(rho2),(r+0.86)*sin(rho2),'**','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho2*(180/pi));
    elseif pmuscle>=1e-2 && pmuscle<1e-1
        text((r+0.86)*cos(rho2),(r+0.86)*sin(rho2),'*','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho2*(180/pi));
    else
        text((r+0.86)*cos(rho2),(r+0.86)*sin(rho2),'n.s.','FontSize',FSize,'HorizontalAlignment','center','rotation',90+rho2*(180/pi));
    end
end
text(0,0,'EMS','HorizontalAlignment','center','FontName','Arial','FontSize',FSize);

% shape descriptor name
text(1.5*cos((1-0.65)*dangle),1.5*sin((1-0.65)*dangle),{'General','Sphericity'},'HorizontalAlignment','center','FontName','Arial','FontSize',FSize);
text(1.4*cos((2-0.6)*dangle),1.4*sin((2-0.8)*dangle),{'Intercept','Sphericity'},'HorizontalAlignment','center','FontName','Arial','FontSize',FSize);
text(1.55*cos((3-0.65)*dangle),1.55*sin((3-0.65)*dangle),{'Hayakawa','Roundness'},'HorizontalAlignment','center','FontName','Arial','FontSize',FSize);
text(1.45*cos((4-0.7)*dangle),1.45*sin((4-0.7)*dangle),{'Elongation','Ratio'},'HorizontalAlignment','center','FontName','Arial','FontSize',FSize);
text(1.5*cos((5-0.5)*dangle),1.5*sin((5-0.5)*dangle),{'Huang','Shape Factor'},'HorizontalAlignment','center','FontName','Arial','FontSize',FSize);

% legend
a=1.8;b=2;c=-0.06;d=0.06;
patch([a,b,b,a],[c+0.2,c+0.2,d+0.2,d+0.2],PharynxColor,'LineWidth',LinW,'FaceAlpha',0.7,'EdgeColor','k');hold on;
text(b+0.08,0.2,'MS: Pharynx','FontSize',FSize,'HorizontalAlignment','left');
patch([a,b,b,a],[c,c,d,d],MuscleColor,'LineWidth',LinW,'FaceAlpha',0.7,'EdgeColor','k');hold on;
text(b+0.08,0,'MS: Muscle','FontSize',FSize,'HorizontalAlignment','left');
patch([a,b,b,a],[c-0.2,c-0.2,d-0.2,d-0.2],IntestineColor,'LineWidth',LinW,'FaceAlpha',0.7,'EdgeColor','k');hold on;
text(b+0.08,-0.2,'E: Intestine','FontSize',FSize,'HorizontalAlignment','left');

axis equal;
ax=1.65;
axis([-ax-0.2,3,-ax,ax]);
axis off;
set(gca,'FontName','Arial','FontSize',FSize);
set(gcf,'units','centimeters','Position',[10 5 25 16]);

