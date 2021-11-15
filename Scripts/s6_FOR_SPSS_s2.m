for CONDI4Corr=2:3;  % Std, Targ, Nov
    IDENTITY.ERP(:,CONDI4Corr-1)=squeeze(mean(MEGA_ERP(:,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3));
end

UNIQUE_SX=unique(IDENTITY.DEMO(~isnan(IDENTITY.DEMO(:,1)),1));
for sxi=1:length(UNIQUE_SX)
    thisguy=UNIQUE_SX(sxi);
    
    FIRST=[]; SECOND=[]; 
    
    FIRST=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==1) ));
    SECOND=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==2) ));

    % for standard models
    if ~isempty(FIRST)  % B/C of bad EEG
        FORSPSS(sxi,1)=IDENTITY.DEMO(FIRST,1);                                          FORSPSS_HDR{1}='subno';
        FORSPSS(sxi,2)=IDENTITY.DEMO(FIRST,find(strcmp('session',IDENTITY_DEMO_HDR)));  FORSPSS_HDR{2}='session';
        FORSPSS(sxi,3)=IDENTITY.DEMO(FIRST,find(strcmp('group',IDENTITY_DEMO_HDR)));    FORSPSS_HDR{3}='group';
        FORSPSS(sxi,4)=IDENTITY.DEMO(FIRST,find(strcmp('SexF1',IDENTITY_DEMO_HDR)));    FORSPSS_HDR{4}='SexF1';
        FORSPSS(sxi,5)=IDENTITY.DEMO(FIRST,find(strcmp('Age',IDENTITY_DEMO_HDR)));      FORSPSS_HDR{5}='Age';
        FORSPSS(sxi,6)=IDENTITY.NP(FIRST,find(strcmp('TOPF',IDENTITY_NP_HDR)));         FORSPSS_HDR{6}='TOPF';
        FORSPSS(sxi,7)=IDENTITY.TBI(FIRST,find(strcmp('Days',IDENTITY_TBI_HDR)));       FORSPSS_HDR{7}='Days';
        
        FORSPSS(sxi,8)=IDENTITY.QUEX(FIRST,find(strcmp('BDI',IDENTITY_QUEX_HDR)));      FORSPSS_HDR{8}='BDI_1';
        FORSPSS(sxi,9)=IDENTITY.QUEX(FIRST,find(strcmp('NSItot',IDENTITY_QUEX_HDR)));   FORSPSS_HDR{9}='NSI_1';
        FORSPSS(sxi,10)=IDENTITY.QUEX(FIRST,find(strcmp('F_Tot_B4',IDENTITY_QUEX_HDR))); FORSPSS_HDR{10}='F_B4_1';
        FORSPSS(sxi,11)=IDENTITY.QUEX(FIRST,find(strcmp('F_Tot',IDENTITY_QUEX_HDR)));   FORSPSS_HDR{11}='F_Tot_1';
        
        FORSPSS(sxi,12)=IDENTITY.ERP(FIRST,1);                                          FORSPSS_HDR{12}='P3b_1';
        FORSPSS(sxi,13)=IDENTITY.ERP(FIRST,2);                                          FORSPSS_HDR{13}='P3a_1';
    else
        FORSPSS(sxi,1)=IDENTITY.DEMO(SECOND,1);                                          
        FORSPSS(sxi,2)=IDENTITY.DEMO(SECOND,find(strcmp('session',IDENTITY_DEMO_HDR)));   
        FORSPSS(sxi,3)=IDENTITY.DEMO(SECOND,find(strcmp('group',IDENTITY_DEMO_HDR)));     
        FORSPSS(sxi,4)=IDENTITY.DEMO(SECOND,find(strcmp('SexF1',IDENTITY_DEMO_HDR)));    
        FORSPSS(sxi,5)=IDENTITY.DEMO(SECOND,find(strcmp('Age',IDENTITY_DEMO_HDR)));      
        FORSPSS(sxi,6)=NaN;
        FORSPSS(sxi,7)=NaN;
        
        FORSPSS(sxi,8)=NaN;
        FORSPSS(sxi,9)=NaN;
        FORSPSS(sxi,10)=NaN;
        FORSPSS(sxi,11)=NaN;
        
        FORSPSS(sxi,12)=NaN;
        FORSPSS(sxi,13)=NaN;       
    end
    
    if ~isempty(SECOND)
        FORSPSS(sxi,14)=IDENTITY.QUEX(SECOND,find(strcmp('BDI',IDENTITY_QUEX_HDR)));      FORSPSS_HDR{14}='BDI_2';
        FORSPSS(sxi,15)=IDENTITY.QUEX(SECOND,find(strcmp('NSItot',IDENTITY_QUEX_HDR)));   FORSPSS_HDR{15}='NSI_2';
        FORSPSS(sxi,16)=IDENTITY.QUEX(SECOND,find(strcmp('F_Tot_B4',IDENTITY_QUEX_HDR))); FORSPSS_HDR{16}='F_B4_2';
        FORSPSS(sxi,17)=IDENTITY.QUEX(SECOND,find(strcmp('F_Tot',IDENTITY_QUEX_HDR)));   FORSPSS_HDR{17}='F_Tot_2';
        FORSPSS(sxi,18)=IDENTITY.ERP(SECOND,1);                                          FORSPSS_HDR{18}='P3b_2';
        FORSPSS(sxi,19)=IDENTITY.ERP(SECOND,2);                                          FORSPSS_HDR{19}='P3a_2';
    else
        FORSPSS(sxi,14)=NaN;
        FORSPSS(sxi,15)=NaN;
        FORSPSS(sxi,16)=NaN;
        FORSPSS(sxi,17)=NaN;
        FORSPSS(sxi,18)=NaN;
        FORSPSS(sxi,19)=NaN;
    end
    
end

% CTL vs. mmTBI
[H,P,CI,STATS]=ttest2(FORSPSS(FORSPSS(:,3)==1,[12,13]),FORSPSS(FORSPSS(:,3)==3,[12,13]))
% subacute vs. mmTBI
[H,P,CI,STATS]=ttest2(FORSPSS(FORSPSS(:,3)==2,[12,13]),FORSPSS(FORSPSS(:,3)==3,[12,13]))
 

 % For Mixed Linear Modeling
 for sxi=1:length(IDENTITY.DEMO)
     
     FORMLM(sxi,1)=IDENTITY.DEMO(sxi,1);                                          FORMLM_HDR{1}='subno';
     FORMLM(sxi,2)=IDENTITY.DEMO(sxi,find(strcmp('session',IDENTITY_DEMO_HDR)));  FORMLM_HDR{2}='session';
     FORMLM(sxi,3)=IDENTITY.DEMO(sxi,find(strcmp('group',IDENTITY_DEMO_HDR)));    FORMLM_HDR{3}='group';
     FORMLM(sxi,4)=IDENTITY.DEMO(sxi,find(strcmp('SexF1',IDENTITY_DEMO_HDR)));    FORMLM_HDR{4}='SexF1';
     FORMLM(sxi,5)=IDENTITY.DEMO(sxi,find(strcmp('Age',IDENTITY_DEMO_HDR)));      FORMLM_HDR{5}='Age';
     FORMLM(sxi,6)=IDENTITY.NP(sxi,find(strcmp('TOPF',IDENTITY_NP_HDR)));         FORMLM_HDR{6}='TOPF';
     FORMLM(sxi,7)=IDENTITY.TBI(sxi,find(strcmp('Days',IDENTITY_TBI_HDR)));       FORMLM_HDR{7}='Days';
     
     FORMLM(sxi,8)=IDENTITY.QUEX(sxi,find(strcmp('BDI',IDENTITY_QUEX_HDR)));      FORMLM_HDR{8}='BDI';
     FORMLM(sxi,9)=IDENTITY.QUEX(sxi,find(strcmp('NSItot',IDENTITY_QUEX_HDR)));   FORMLM_HDR{9}='NSI';
     FORMLM(sxi,10)=IDENTITY.QUEX(sxi,find(strcmp('F_Tot_B4',IDENTITY_QUEX_HDR))); FORMLM_HDR{10}='F_B4';
     FORMLM(sxi,11)=IDENTITY.QUEX(sxi,find(strcmp('F_Tot',IDENTITY_QUEX_HDR)));   FORMLM_HDR{11}='F_Tot';
     
     FORMLM(sxi,12)=IDENTITY.ERP(sxi,1);                                          FORMLM_HDR{12}='P3b';
     FORMLM(sxi,13)=IDENTITY.ERP(sxi,2);                                          FORMLM_HDR{13}='P3a';
     
 end

FORMLM_HDR=FORMLM_HDR';
FORSPSS_HDR=FORSPSS_HDR';

