clear all;clc;close all;

load('EMSGeneList.mat');
PharynxGene=EMSGeneList{1,1};
MuscleGene=EMSGeneList{1,2};
PhaMusGene=EMSGeneList{1,3};
IntestineGene=EMSGeneList{1,4};

% fontsize according to rank
gl=16;
FontList=mapminmax([gl:-1:1],5,16);
PharynxGene(1:gl,2)=num2cell(flipud(FontList));
PharynxGene(gl+1:end,2)=num2cell(repmat(5,size(PharynxGene,1)-gl,1));
gl=20;
FontList=mapminmax([gl:-1:1],5.5,16);
PhaMusGene(1:gl,2)=num2cell(flipud(FontList));
PhaMusGene(gl+1:end,2)=num2cell(repmat(5.5,size(PhaMusGene,1)-gl,1));
gl=22;
FontList=mapminmax([gl:-1:1],5.5,16);
MuscleGene(1:gl,2)=num2cell(flipud(FontList));
MuscleGene(gl+1:end,2)=num2cell(repmat(5.5,size(MuscleGene,1)-gl,1));
gl=10;
FontList=mapminmax([gl:-1:1],5,14);
IntestineGene(1:gl,2)=num2cell(flipud(FontList));
IntestineGene(gl+1:end,2)=num2cell(repmat(5,size(IntestineGene,1)-gl,1));

% Venn
LinW=2;
PharynxLine=[132,62,160]./255;
PharynxPat=[158,119,221]./255;
MuscleLine=[231,121,184]./255;
MusclePat=[235,143,182]./255;
IntestineLine=[235,194,61]./255;
IntestinePat=[250,199,86]./255;
% radius
rp=sqrt(length(PharynxGene)+length(PhaMusGene));
rm=sqrt(length(MuscleGene)+length(PhaMusGene));
ri=sqrt(length(IntestineGene));
% center coordinate
xp=0;yp=0;
xm=13;ym=0;
xi=42;yi=0;

% angel [0,2*pi]
theta=0:pi/40:2*pi;
% Pharynx gene
xxp=xp+rp*cos(theta);
yyp=yp+rp*sin(theta);
fill(xxp,yyp,PharynxPat,'FaceAlpha',0.3,'EdgeColor','none');hold on;
% Muscle gene
xxm=xm+rm*cos(theta);
yym=ym+rm*sin(theta);
fill(xxm,yym,MusclePat,'FaceAlpha',0.3,'EdgeColor','none');hold on;
% Intestine gene
xxi=xi+ri*cos(theta);
yyi=yi+ri*sin(theta);
fill(xxi,yyi,IntestinePat,'FaceAlpha',0.3,'EdgeColor','none');hold on;

% pharynx: dash line
anglepha=acos( (rp^2+(xm-xp)^2-rm^2)./(2*rp*(xm-xp)) );
xphadot=xp+rp*cos([anglepha:pi/100:2*pi-anglepha]);
yphadot=yp+rp*sin([anglepha:pi/100:2*pi-anglepha]);
plot(xphadot,yphadot,'--','Color',PharynxLine,'LineWidth',LinW);hold on;
% pharynx: full line
xphafull=xp+rp*cos([-anglepha:pi/100:anglepha]);
yphafull=yp+rp*sin([-anglepha:pi/100:anglepha]);
plot(xphafull,yphafull,'-','Color',PharynxLine,'LineWidth',LinW);hold on;
% muscle: dash line
anglemus=acos( (rm^2+(xm-xp)^2-rp^2)./(2*rm*(xm-xp)) );
xmusdot=xm+rm*cos([-pi+anglemus:pi/100:pi-anglemus]);
ymusdot=ym+rm*sin([-pi+anglemus:pi/100:pi-anglemus]);
plot(xmusdot,ymusdot,'--','Color',MuscleLine,'LineWidth',LinW);hold on;
% muscle: full line
xmusfull=xm+rm*cos([pi-anglemus:pi/100:pi+anglemus]);
ymusfull=ym+rm*sin([pi-anglemus:pi/100:pi+anglemus]);
plot(xmusfull,ymusfull,'-','Color',MuscleLine,'LineWidth',LinW);hold on;
% intestine: full line
plot(xxi,yyi,'-','Color',IntestineLine,'LineWidth',LinW);hold on;

axis equal off;
set(gcf,'units','centimeters','Position',[10,5,22,16]);

% write gene
bas=40;randnum=20000;

% --- Helper to test overlap of two boxes a and b ---
doOverlap = @(a,b) ~( (a(3) < b(1)) || (b(3) < a(1)) ...
    || (a(4) < b(2)) || (b(4) < a(2)) );

% pharynx,40 genes
h=text(-10,0,num2str(size(PharynxGene,1)),'FontName','Arial','FontSize',18,'BackgroundColor','w','FontWeight','bold');
ext=h.Extent;pad=ext(4)/bas;
t_ext=[ext(1)-10*pad,ext(2)-10*pad,ext(1)+ext(3)+20*pad,ext(2)+ext(4)+20*pad];
existExt=t_ext;
for i=1:40
    pr=rand(randnum,2);anglemus=acos( (rm^2+(xm-xp)^2-rp^2)./(2*rm*(xm-xp)) );
    % sprinkle pharynx circle for 1000 times
    for j=1:size(pr,1)
        rpj=rp*pr(j,1);thetaj=2*pi*pr(j,2);
        x=xp+rpj*cos(thetaj);y=yp+rpj*sin(thetaj);
        h=text(x,y,['\it',PharynxGene{i,1}],'FontName','Arial','FontSize',PharynxGene{i,2});
        drawnow;
        ext=h.Extent;pad=ext(4)/bas;
        % text boundary
        t_ext=[ext(1)-2*pad,ext(2)+pad,ext(1)+ext(3)+4*pad,ext(2)+ext(4)-2*pad];
        corners1=[t_ext(1)-xp,t_ext(2)-yp;t_ext(3)-xp,t_ext(2)-yp;t_ext(1)-xp,t_ext(4)-yp;t_ext(3)-xp,t_ext(4)-yp];
        corners2=[t_ext(1)-xm,t_ext(2)-ym;t_ext(3)-xm,t_ext(2)-ym;t_ext(1)-xm,t_ext(4)-ym;t_ext(3)-xm,t_ext(4)-ym];
        % no overlap
        if all(sum(corners1.^2,2)<=(rp-0.1)^2) && all(sum(corners2.^2,2)>(rm+0.1)^2)
            % check no overlap with previous boxes
            if isempty(existExt) || ~any( arrayfun(@(i) doOverlap(existExt(i,:), t_ext), 1:size(existExt,1)) )
                % success!
                existExt(end+1,:) = t_ext;
                break;
            end
        end
        % if we get here, fail: delete and retry
        delete(h);
    end
    if j==size(pr,1)
        warning('Could not place "%s" after %d tries.', PharynxGene{i,1}, j);
    end
end

% muscle,60 genes
h=text(19,0,num2str(size(MuscleGene,1)),'FontName','Arial','FontSize',18,'BackgroundColor','w','FontWeight','bold');
ext=h.Extent;pad=ext(4)/bas;
t_ext=[ext(1)-10*pad,ext(2)-10*pad,ext(1)+ext(3)+20*pad,ext(2)+ext(4)+20*pad];
existExt=t_ext;
for i=1:100
    pm=rand(randnum,2);
    % sprinkle muscle circle for 1000 times
    for j=1:size(pm,1)
        rmj=rm*pm(j,1);thetaj=2*pi*pm(j,2);
        x=xm+rmj*cos(thetaj);y=ym+rmj*sin(thetaj);
        h=text(x,y,['\it',MuscleGene{i,1}],'FontName','Arial','FontSize',MuscleGene{i,2});
        drawnow;
        ext=h.Extent;pad=ext(4)/bas;
        % text boundary
        t_ext=[ext(1)-2*pad,ext(2)+pad,ext(1)+ext(3)+4*pad,ext(2)+ext(4)-2*pad];
        corners1=[t_ext(1)-xp,t_ext(2)-yp;t_ext(3)-xp,t_ext(2)-yp;t_ext(1)-xp,t_ext(4)-yp;t_ext(3)-xp,t_ext(4)-yp];
        corners2=[t_ext(1)-xm,t_ext(2)-ym;t_ext(3)-xm,t_ext(2)-ym;t_ext(1)-xm,t_ext(4)-ym;t_ext(3)-xm,t_ext(4)-ym];
        % no overlap
        if all(sum(corners1.^2,2)>(rp+0.1)^2) && all(sum(corners2.^2,2)<=(rm-0.1)^2)
            % check no overlap with previous boxes
            if isempty(existExt) || ~any( arrayfun(@(i) doOverlap(existExt(i,:), t_ext), 1:size(existExt,1)) )
                % success!
                existExt(end+1,:) = t_ext;
                break;
            end
        end
        % if we get here, fail: delete and retry
        delete(h);
    end
    if j==size(pm,1)
        warning('Could not place "%s" after %d tries.', MuscleGene{i,1}, j);
    end
end

% pharynx & Muscle,40 genes
h=text(4,0,num2str(size(PhaMusGene,1)),'FontName','Arial','FontSize',18,'BackgroundColor','w','FontWeight','bold');
ext=h.Extent;pad=ext(4)/bas;
t_ext=[ext(1)-10*pad,ext(2)-10*pad,ext(1)+ext(3)+20*pad,ext(2)+ext(4)+20*pad];
existExt=t_ext;
for i=1:100
    pr=rand(randnum,2);
    % sprinkle pharynx and muscle overlapping circle for 1000 times
    for j=1:size(pr,1)
        rpj=rp*pr(j,1);thetaj=2*pi*pr(j,2);
        x=xp+rpj*cos(thetaj);y=yp+rpj*sin(thetaj);
        h=text(x,y,['\it',PhaMusGene{i,1}],'FontName','Arial','FontSize',PhaMusGene{i,2});
        drawnow;
        ext=h.Extent;pad=ext(4)/bas;
        % text boundary
        t_ext=[ext(1)-2*pad,ext(2)+pad,ext(1)+ext(3)+4*pad,ext(2)+ext(4)-2*pad];
        corners1=[t_ext(1)-xp,t_ext(2)-yp;t_ext(3)-xp,t_ext(2)-yp;t_ext(1)-xp,t_ext(4)-yp;t_ext(3)-xp,t_ext(4)-yp];
        corners2=[t_ext(1)-xm,t_ext(2)-ym;t_ext(3)-xm,t_ext(2)-ym;t_ext(1)-xm,t_ext(4)-ym;t_ext(3)-xm,t_ext(4)-ym];
        % no overlap
        if all(sum(corners1.^2,2)<=(rp-0.1)^2) && all(sum(corners2.^2,2)<=(rm-0.1)^2)
            % check no overlap with previous boxes
            if isempty(existExt) || ~any( arrayfun(@(i) doOverlap(existExt(i,:), t_ext), 1:size(existExt,1)) )
                % success!
                existExt(end+1,:) = t_ext;
                break;
            end
        end
        % if we get here, fail: delete and retry
        delete(h);
    end
    if j==size(pm,1)
        warning('Could not place "%s" after %d tries.', PhaMusGene{i,1}, j);
    end
end

% intestine,50 genes
h=text(41,0,num2str(size(IntestineGene,1)),'FontName','Arial','FontSize',18,'BackgroundColor','w','FontWeight','bold');
ext=h.Extent;pad=ext(4)/bas;
t_ext=[ext(1)-10*pad,ext(2)-10*pad,ext(1)+ext(3)+20*pad,ext(2)+ext(4)+20*pad];
existExt=t_ext;
for i=1:50
    pit=rand(randnum,2);
    % sprinkle intestine overlapping circle for 1000 times
    for j=1:size(pm,1)
        rij=ri*pit(j,1);thetaj=2*pi*pit(j,2);
        x=xi+rij*cos(thetaj);y=yi+rij*sin(thetaj);
        h=text(x,y,['\it',IntestineGene{i,1}],'FontName','Arial','FontSize',IntestineGene{i,2});
        drawnow;
        ext=h.Extent;pad=ext(4)/bas;
        % text boundary
        t_ext=[ext(1)-2*pad,ext(2)+pad,ext(1)+ext(3)+4*pad,ext(2)+ext(4)-2*pad];
        corners=[t_ext(1)-xi,t_ext(2)-yi;t_ext(3)-xi,t_ext(2)-yi;t_ext(1)-xi,t_ext(4)-yi;t_ext(3)-xi,t_ext(4)-yi];
        % no overlap
        if all(sum(corners.^2,2)<=ri^2)
            % check no overlap with previous boxes
            if isempty(existExt) || ~any( arrayfun(@(i) doOverlap(existExt(i,:), t_ext), 1:size(existExt,1)) )
                % success!
                existExt(end+1,:) = t_ext;
                break;
            end
        end
        % if we get here, fail: delete and retry
        delete(h);
    end
    if j==size(pit,1)
        warning('Could not place "%s" after %d tries.', IntestineGene{i,1}, j);
    end
end

% save image
saveas(gcf,['EMS_Venn3.svg'],'svg');
print(gcf,['EMS_Venn3'],'-dpng','-r600');
