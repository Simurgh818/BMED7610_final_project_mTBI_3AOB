
cd(datadir);

filz=dir(['*_3AOB_TFandERPs_V.mat']);
Nsubjs=length(filz);

% Load BigAgg_Data
load('C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\Scripts\BigAgg_Data.mat');

% Preallocate
IDENTITY.DEMO=NaN(Nsubjs,7);
IDENTITY.TBI=NaN(Nsubjs,9);
IDENTITY.NP=NaN(Nsubjs,8);
IDENTITY.QUEX=NaN(Nsubjs,25);
% ^^^^^^^^^^
MEGA_PWR=NaN(Nsubjs,3,50,751,3);
MEGA_PHS=NaN(Nsubjs,3,50,751,3);
MEGA_SYNCH=NaN(Nsubjs,3,50,751,3);
MEGA_SYNCH_TOPO=NaN(Nsubjs,60,751,3);
MEGA_POWER_TOPO=NaN(Nsubjs,60,751,3);
MEGA_CORREL_TOPO=NaN(Nsubjs,60,751,3);
MEGA_ERP=NaN(Nsubjs,60,751,3);
MEGA_TRL_ct=NaN(Nsubjs,3);

for si=1:Nsubjs
    
    subno = str2double(filz(si).name(1:end-23)) ; % B/C some Quinn ones have 5 digit IDs, some 4
    session = str2double(filz(si).name(end-21)) ;
    
    if subno<3500 % Cavanagh
        if mod(subno,2)==1, group=1;     % ODD - Ctl
        elseif mod(subno,2)==0,group=2;  % EVEN - mTBI
        end
    else % Quinn
        group=3;    % Chronic TBI  (cTBI)
    end
    
    IDENTITY_DEMO_HDR{1}={'subno'}; IDENTITY_DEMO_HDR{2}='session'; IDENTITY_DEMO_HDR{3}='group'; 
    IDENTITY_TBI_HDR{1}={'subno'};  IDENTITY_TBI_HDR{2}='session';  IDENTITY_TBI_HDR{3}='group'; 
    IDENTITY_NP_HDR{1}={'subno'};   IDENTITY_NP_HDR{2}='session';   IDENTITY_NP_HDR{3}='group'; 
    IDENTITY_QUEX_HDR{1}={'subno'}; IDENTITY_QUEX_HDR{2}='session'; IDENTITY_QUEX_HDR{3}='group'; 
    
    IDENTITY.DEMO(si,1:3)=[subno,session,group];
    IDENTITY.TBI(si,1:3)=[subno,session,group];
    IDENTITY.NP(si,1:3)=[subno,session,group];
    IDENTITY.QUEX(si,1:3)=[subno,session,group];
    
    % --------------- QUEX
    if group<3  % Cavanagh data
        if any(DEMO.ID(:,1)==subno)
            bigagg_idx=find(DEMO.ID(:,1)==subno);
                IDENTITY.DEMO(si,4)=DEMO.URSI(bigagg_idx,1); IDENTITY_DEMO_HDR{4}='URSI';
                IDENTITY.DEMO(si,5)=DEMO.Sex_F1(bigagg_idx); IDENTITY_DEMO_HDR{5}='SexF1';
                IDENTITY.DEMO(si,6)=DEMO.Age(bigagg_idx);    IDENTITY_DEMO_HDR{6}='Age';
                IDENTITY.DEMO(si,7)=DEMO.SES(bigagg_idx);    IDENTITY_DEMO_HDR{7}='YrsEd';  % That's what this actually is.
            if session==1
                IDENTITY.TBI(si,4)=TBIfields.Glasgow(bigagg_idx);          IDENTITY_TBI_HDR{4}='GCS';
                IDENTITY.TBI(si,5)=TBIfields.LOC(bigagg_idx);              IDENTITY_TBI_HDR{5}='LOC';
                IDENTITY.TBI(si,6)=TBIfields.LOCtime(bigagg_idx);          IDENTITY_TBI_HDR{6}='LOCtime';
                IDENTITY.TBI(si,7)=TBIfields.LOM(bigagg_idx);              IDENTITY_TBI_HDR{7}='LOM';
                IDENTITY.TBI(si,8)=TBIfields.DaysSinceInjury(bigagg_idx);  IDENTITY_TBI_HDR{8}='Days';
                % ---------------------
                IDENTITY.NP(si,4)=NP.TOPF_Score(bigagg_idx);               IDENTITY_NP_HDR{4}='TOPF';
                IDENTITY.NP(si,5)=NP.Coding(bigagg_idx);                   IDENTITY_NP_HDR{5}='Coding';
                IDENTITY.NP(si,6)=NP.SPAN.Tot(bigagg_idx);                 IDENTITY_NP_HDR{6}='Span';
                IDENTITY.NP(si,7)=mean([NP.HVLT.T1(bigagg_idx),NP.HVLT.T2(bigagg_idx),NP.HVLT.T3(bigagg_idx)]'); IDENTITY_NP_HDR{7}='HVLT13';
                IDENTITY.NP(si,8)=NP.HVLT.DelayRecall(bigagg_idx);         IDENTITY_NP_HDR{8}='HVLTdelay';
            end
            IDENTITY.QUEX(si,4)=QUEX.BDI(bigagg_idx,session);           IDENTITY_QUEX_HDR{4}='BDI';
            IDENTITY.QUEX(si,5)=QUEX.NSI.tot(bigagg_idx,session);       IDENTITY_QUEX_HDR{5}='NSItot';
            IDENTITY.QUEX(si,6)=QUEX.NSI.somatic(bigagg_idx,session);   IDENTITY_QUEX_HDR{6}='NSIsom';
            IDENTITY.QUEX(si,7)=QUEX.NSI.cog(bigagg_idx,session);       IDENTITY_QUEX_HDR{7}='NSIcog';
            IDENTITY.QUEX(si,8)=QUEX.NSI.emo(bigagg_idx,session);       IDENTITY_QUEX_HDR{8}='NSIemo';
            IDENTITY.QUEX(si,9)=QUEX.FRSBE.Tot_B4(bigagg_idx,session); IDENTITY_QUEX_HDR{9}='F_Tot_B4';   
            IDENTITY.QUEX(si,10)=QUEX.FRSBE.Tot_Now(bigagg_idx,session); IDENTITY_QUEX_HDR{10}='F_Tot';
            IDENTITY.QUEX(si,11)=EX.EX(bigagg_idx,session);   IDENTITY_QUEX_HDR{11}='EX_EX';
            IDENTITY.QUEX(si,12)=EX.CC(bigagg_idx,session);   IDENTITY_QUEX_HDR{12}='EX_CC';
            IDENTITY.QUEX(si,13)=EX.FL(bigagg_idx,session);   IDENTITY_QUEX_HDR{13}='EX_FL';
            IDENTITY.QUEX(si,14)=EX.WM(bigagg_idx,session);   IDENTITY_QUEX_HDR{14}='EX_WM';
        end
    elseif group==3
        if any(Q_DEMO.URSI==subno)
            bigagg_idx=find(Q_DEMO.URSI==subno);
            %^^^^^^^^^^^
            IDENTITY.DEMO(si,4)=Q_DEMO.URSI(bigagg_idx,1);
            IDENTITY.DEMO(si,5)=Q_DEMO.Sex_F1(bigagg_idx);
            IDENTITY.DEMO(si,6)=Q_DEMO.Age(bigagg_idx);
            IDENTITY.DEMO(si,7)=Q_DEMO.SES(bigagg_idx);
            % ---------------------
            IDENTITY.TBI(si,6)=Q_TBIfields.LOCdurMINS(bigagg_idx);
            IDENTITY.TBI(si,9)=Q_TBIfields.YearsSinceInjury(bigagg_idx);  IDENTITY_TBI_HDR{9}='Years';
            % ---------------------
            IDENTITY.NP(si,4)=Q_NP.TOPF_Score(bigagg_idx);
            IDENTITY.NP(si,5)=Q_NP.Coding(bigagg_idx);
            IDENTITY.NP(si,6)=Q_NP.SPAN.Tot(bigagg_idx);
            IDENTITY.NP(si,7)=mean([Q_NP.HVLT.T1(bigagg_idx),Q_NP.HVLT.T2(bigagg_idx),Q_NP.HVLT.T3(bigagg_idx)]');
            IDENTITY.NP(si,8)=Q_NP.HVLT.DelayRecall(bigagg_idx);
            % ---------------------
            IDENTITY.QUEX(si,4)=Q_QUEX.BDI(bigagg_idx);            
            IDENTITY.QUEX(si,5)=Q_QUEX.NSI.tot(bigagg_idx);        
            IDENTITY.QUEX(si,6)=Q_QUEX.NSI.somatic(bigagg_idx);   
            IDENTITY.QUEX(si,7)=Q_QUEX.NSI.cog(bigagg_idx);        
            IDENTITY.QUEX(si,8)=Q_QUEX.NSI.emo(bigagg_idx);        
            IDENTITY.QUEX(si,9)=Q_QUEX.FRSBE.RAW_Tot_Now(bigagg_idx);     
            IDENTITY.QUEX(si,10)=Q_QUEX.FRSBE.Tot_Now(bigagg_idx);      
        end
    end
    clear bigagg_idx
    
    % EEG
    load([filz(si).name(1:end-5),'V.mat'],'ERP','TRL_ct');  
    MEGA_ERP(si,:,:,:)=ERP;
    MEGA_TRL_ct(si,:)=TRL_ct;

    clear ERP ISPC ITPC POWER VECTOR RT SYNCH_TOPO TRL_ct subno session group BEH ACC RT POWER_TOPO CORREL_TOPO;
    
end

clear DEMO QUEX NP EX TBIfields Q_* 

cd(homedir);
