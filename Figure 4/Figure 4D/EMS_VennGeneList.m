% EMS Venn Diagram
clear all; clc; close all;

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

PharynxGene={};
MuscleGene={};
PhaMusGene={};
IntestineGene={};
for i=1:size(EMSAll,2)
    genename=EMSAll{1,i};
    % p value and asymmetry of all and partial cells
    PsigPar1=EMSPartial1{5,i};
    AsymPar1=EMSPartial1{6,i};
    PsigPar2=EMSPartial2{5,i};
    AsymPar2=EMSPartial2{6,i};
    PsigAll=EMSAll{5,i};
    AsymAll=EMSAll{6,i};
    % pharynx differentiation
    if PsigPar1<0.1
        % single strengthen
        if (PsigPar1<PsigAll) || (abs(AsymPar1)>=abs(AsymAll))
            if AsymPar1>0
                % Pharynx marker
                PharynxGene=[PharynxGene;{genename,floor(-log10(PsigPar1)),abs(AsymPar1)}];
            else
                IntestineGene=[IntestineGene;{genename,floor(-log10(PsigPar1)),abs(AsymPar1)}];
            end
        end
    end
    % muscle differentiation
    if PsigPar2<0.1
        % single strengthen
        if (PsigPar2<PsigAll) || (abs(AsymPar2)>=abs(AsymAll))
            if AsymPar2>0
                % Muscle marker
                MuscleGene=[MuscleGene;{genename,floor(-log10(PsigPar2)),abs(AsymPar2)}];
            else
                IntestineGene=[IntestineGene;{genename,floor(-log10(PsigPar2)),abs(AsymPar2)}];
            end
        end
    end
end
% sort
[sortPha,IP]=sortrows(cell2mat(PharynxGene(:,[2,3])),[1,2],'descend');
[sortMus,IM]=sortrows(cell2mat(MuscleGene(:,[2,3])),[1,2],'descend');
[sortIns,II]=sortrows(cell2mat(IntestineGene(:,[2,3])),[1,2],'descend');
SPharynxGene=PharynxGene(IP,1);
SMuscleGene=MuscleGene(IM,1);
SIntestineGene=IntestineGene(II,1);

% shared genes
PhaMusGene=intersect(SPharynxGene,SMuscleGene,'stable');
IntestineGene=unique(SIntestineGene,'stable');
PharynxGene=setdiff(SPharynxGene,PhaMusGene,'stable');
MuscleGene=setdiff(SMuscleGene,PhaMusGene,'stable');

EMSGeneList={PharynxGene,MuscleGene,PhaMusGene,IntestineGene};
save('EMSGeneList.mat','EMSGeneList','-v7.3');

