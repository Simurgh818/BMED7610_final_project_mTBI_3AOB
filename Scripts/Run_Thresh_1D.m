function [Corrected_P] = Run_Thresh_1D(TEMP1,TEMP2,site,ttesttype)

SHUFFLES=5000;
for shuffi=1:SHUFFLES
    
    tempA=squeeze(mean(TEMP1(:,site,251:1001),2));
    tempB=squeeze(mean(TEMP2(:,site,251:1001),2));
    tempAB=[tempA;tempB];
    idx=shuffle([ones(1,25),zeros(1,25)]);
    A=tempAB(idx==1,:);
    B=tempAB(idx==0,:);
    
    if strmatch(ttesttype,'between')
         [H,P,CI,STATS]=ttest2(A,B);
    elseif strmatch(ttesttype,'within')
         [H,P,CI,STATS]=ttest(A,B);
    end
    P(P<=.05)=NaN; P(P>.05)=0; P(isnan(P))=1;
    
    P=squeeze(P);
    l=bwlabel(P);
    if max(l)>0
        for ei=1:max(l)
            idxs = find(l == ei);
            tempthresh(ei) = sum(abs(STATS.tstat(idxs)));
        end
        THRESH(shuffi) = max(tempthresh);
        clear idxs tempthresh ;
    else
        THRESH(shuffi) = 0;
    end
    clear H CI P STATS temp* A B idx  l dims lmax;
    
end
THRESH=sort(THRESH);
ThisThreshold=THRESH(end-SHUFFLES*.05);


% NOW Run 1D size of effects
if strmatch(ttesttype,'between')
    [H,P,CI,STATS]=ttest2(squeeze(mean(TEMP1(:,site,251:1001),2)),squeeze(mean(TEMP2(:,site,251:1001),2)));
elseif strmatch(ttesttype,'within')
    [H,P,CI,STATS]=ttest(squeeze(mean(TEMP1(:,site,251:1001),2)),squeeze(mean(TEMP2(:,site,251:1001),2)));
end
P(P<=.05)=NaN; P(P>.05)=0; P(isnan(P))=1;

P=squeeze(P);
l=bwlabel(P);
Corrected_P=NaN*ones(1,751);
if max(l)>0
    for ei=1:max(l)
        idxs = find(l == ei);
        if sum(abs(STATS.tstat(idxs))) > ThisThreshold
            Corrected_P(idxs) = 1;
        end
    end
end
clear H CI P STATS temp* A B idx  l dims lmax idxs;

clear THRESH ThisThreshold

