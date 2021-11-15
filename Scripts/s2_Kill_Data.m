% Kill Quinn S2
for si=1:length(IDENTITY.DEMO)
    if IDENTITY.DEMO(si,3)==3 && IDENTITY.DEMO(si,2)==2
        IDENTITY.DEMO(si,:)=NaN; IDENTITY.TBI(si,:)=NaN; IDENTITY.NP(si,:)=NaN; IDENTITY.QUEX(si,:)=NaN;
% %         MEGA_PWR(si,:,:,:,:)=NaN;
% %         MEGA_PHS(si,:,:,:,:)=NaN;
% %         MEGA_SYNCH(si,:,:,:,:)=NaN;
% %         MEGA_SYNCH_TOPO(si,:,:,:)=NaN;
        MEGA_ERP(si,:,:,:)=NaN;
    end
end

% Kill malingering Quinn patient  # 43047 (is already excluded... but this will make sure!)
if any(IDENTITY.DEMO(:,1)==43047); BOOM; end

% Kill any Quinn patients who were Cavanagh patients
% F48 3032(Cav URSI: 30454; Quinn URSI: 35957) & F22 3004(Cav URSI: 69117; Quinn URSI: 48880)
badidx=find(IDENTITY.DEMO(:,1)==35957);
IDENTITY.DEMO(badidx,:)=NaN; IDENTITY.TBI(badidx,:)=NaN; IDENTITY.NP(badidx,:)=NaN; IDENTITY.QUEX(badidx,:)=NaN; clear badidx; 
% The other was 3004, who is killed below due to no LOC

% Kill 2 mTBI with out LOC          3004, 3056 == no LOC
badidx=find(IDENTITY.DEMO(:,1)==3004);
IDENTITY.DEMO(badidx,:)=NaN; IDENTITY.TBI(badidx,:)=NaN; IDENTITY.NP(badidx,:)=NaN; IDENTITY.QUEX(badidx,:)=NaN; clear badidx; 
badidx=find(IDENTITY.DEMO(:,1)==3056);
IDENTITY.DEMO(badidx,:)=NaN; IDENTITY.TBI(badidx,:)=NaN; IDENTITY.NP(badidx,:)=NaN; IDENTITY.QUEX(badidx,:)=NaN; clear badidx; 

% Kill people with pre-existing head injuries
badidx=find(IDENTITY.DEMO(:,1)==3024);
IDENTITY.DEMO(badidx,:)=NaN; IDENTITY.TBI(badidx,:)=NaN; IDENTITY.NP(badidx,:)=NaN; IDENTITY.QUEX(badidx,:)=NaN; clear badidx; 

% Kill *sessions* if they had an intervening head injury
for si=1:length(IDENTITY.DEMO)
    if IDENTITY.DEMO(si,1)==3034 && IDENTITY.DEMO(si,2)==3
        badidx=si;
    end
end
% IDENTITY.DEMO(badidx,:)=NaN; IDENTITY.TBI(badidx,:)=NaN; IDENTITY.NP(badidx,:)=NaN; IDENTITY.QUEX(badidx,:)=NaN; clear badidx; 

% Kill people with TOMM score<45
badidx=find(IDENTITY.DEMO(:,1)==14000);
IDENTITY.DEMO(badidx,:)=NaN; IDENTITY.TBI(badidx,:)=NaN; IDENTITY.NP(badidx,:)=NaN; IDENTITY.QUEX(badidx,:)=NaN; clear badidx; 

% Assessment > 2 weeks
badidx=find(IDENTITY.TBI(:,8)>14);
IDENTITY.DEMO(badidx,:)=NaN; IDENTITY.TBI(badidx,:)=NaN; IDENTITY.NP(badidx,:)=NaN; IDENTITY.QUEX(badidx,:)=NaN; clear badidx; 








