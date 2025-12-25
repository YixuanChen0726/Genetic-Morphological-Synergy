% ABplpapp WT contact map
clear all;clc;close all;

CellName='ABplpapp';
AllContactCell=[];
% WT
for SampleNum=1:8
    RawContact=readcell(['..\..\bin\CMap Data\Dataset E\WT_Sample',num2str(SampleNum),'\WT_Sample',num2str(SampleNum),'_Stat.csv']);
    Row1=find(strcmp(RawContact(:,1),CellName));
    Row2=find(strcmp(RawContact(:,2),CellName));
    CellList=[RawContact(Row1,2);RawContact(Row2,1)];
    AllContactCell=[AllContactCell;CellList];
end
% lag-1
for SampleNum=1:2
    RawContact=readcell(['..\..\bin\CMap Data\Dataset E\MT_lag-1_Sample',num2str(SampleNum),'\MT_lag-1_Sample',num2str(SampleNum),'_Stat.csv']);
    Row1=find(strcmp(RawContact(:,1),CellName));
    Row2=find(strcmp(RawContact(:,2),CellName));
    CellList=[RawContact(Row1,2);RawContact(Row2,1)];
    AllContactCell=[AllContactCell;CellList];
end
% CShaper
for SampleNum=4:20
    RawContact=readcell(['..\..\bin\CShaper Data\Sample',num2str(SampleNum,'%02d'),'\Sample',num2str(SampleNum,'%02d'),'_Stat.csv']);
    Row1=find(strcmp(RawContact(:,1),CellName));
    Row2=find(strcmp(RawContact(:,2),CellName));
    CellList=[RawContact(Row1,2);RawContact(Row2,1)];
    AllContactCell=[AllContactCell;CellList];
end
AllContactCell=unique(AllContactCell);

% table:sample/contact area
AllContact=cell(8+2+17+1,length(AllContactCell)+1);
AllContact(1,2:end)=AllContactCell';
% WT
for SampleNum=1:8
    AllContact{SampleNum+1,1}=['WT: Sample',num2str(SampleNum)];
    RawContact=readcell(['..\..\bin\CMap Data\Dataset E\WT_Sample',num2str(SampleNum),'\WT_Sample',num2str(SampleNum),'_Stat.csv']);
    for CellNum=2:size(AllContact,2)
        ContactCell=AllContact{1,CellNum};
        % contact relation
        Row1=intersect(find(strcmp(RawContact(:,1),CellName)),find(strcmp(RawContact(:,2),ContactCell)));
        Row2=intersect(find(strcmp(RawContact(:,2),CellName)),find(strcmp(RawContact(:,1),ContactCell)));
        Row=[Row1;Row2];
        if ~isempty(Row)
            ContactData=RawContact(Row,3:end);
            Col=find(cellfun(@(x) ~ismissing(x) && x>0,ContactData))+2;
            if ~isempty(Col)
                AllContact{SampleNum+1,CellNum}=mean(cell2mat(RawContact(Row,Col)));
            else
                AllContact{SampleNum+1,CellNum}=0;
            end
        else
            AllContact{SampleNum+1,CellNum}=0;
        end
    end
end

% CShaper
for SampleNum=4:20
    AllContact{SampleNum+6,1}=['CShaper: Sample',num2str(SampleNum)];
    RawContact=readcell(['..\..\bin\CShaper Data\Sample',num2str(SampleNum,'%02d'),'\Sample',num2str(SampleNum,'%02d'),'_Stat.csv']);
    for CellNum=2:size(AllContact,2)
        ContactCell=AllContact{1,CellNum};
        Row1=intersect(find(strcmp(RawContact(:,1),CellName)),find(strcmp(RawContact(:,2),ContactCell)));
        Row2=intersect(find(strcmp(RawContact(:,2),CellName)),find(strcmp(RawContact(:,1),ContactCell)));
        Row=[Row1;Row2];
        if ~isempty(Row)
            ContactData=RawContact(Row,3:end);
            Col=find(cellfun(@(x) ~ismissing(x) && x>0,ContactData))+2;
            if ~isempty(Col)
                AllContact{SampleNum+6,CellNum}=mean(cell2mat(RawContact(Row,Col)));
            else
                AllContact{SampleNum+6,CellNum}=0;
            end
        else
            AllContact{SampleNum+6,CellNum}=0;
        end
    end
end

% lag-1
for SampleNum=1:2
    AllContact{SampleNum+26,1}=['lag-1^-: Sample',num2str(SampleNum)];
    RawContact=readcell(['..\..\bin\CMap Data\Dataset E\MT_lag-1_Sample',num2str(SampleNum),'\MT_lag-1_Sample',num2str(SampleNum),'_Stat.csv']);
    for CellNum=2:size(AllContact,2)
        ContactCell=AllContact{1,CellNum};
        % 找到接触关系所在行
        Row1=intersect(find(strcmp(RawContact(:,1),CellName)),find(strcmp(RawContact(:,2),ContactCell)));
        Row2=intersect(find(strcmp(RawContact(:,2),CellName)),find(strcmp(RawContact(:,1),ContactCell)));
        Row=[Row1;Row2];
        if ~isempty(Row)
            ContactData=RawContact(Row,3:end);
            Col=find(cellfun(@(x) ~ismissing(x) && x>0,ContactData))+2;
            if ~isempty(Col)
                AllContact{SampleNum+26,CellNum}=mean(cell2mat(RawContact(Row,Col)));
            else
                AllContact{SampleNum+26,CellNum}=0;
            end
        else
            AllContact{SampleNum+26,CellNum}=0;
        end
    end
end

save([CellName,'_All_ContactRelation.mat'],'AllContact','-v7.3');
