% #######################################################################################################

% Count
for groupi=1:2
    for time=1:3
        Sx=logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==groupi)  );
        TABLE1_COUNT(groupi,time,:)=[sum(Sx),nansum(IDENTITY.DEMO(Sx,5))];
    end
end

% Other neat stuff
for groupi=1:3
    time=1;
    Sx=logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==groupi)  );
    TABLE1_VARS(1,groupi,:)=[mean(IDENTITY.DEMO(Sx,6)),std(IDENTITY.DEMO(Sx,6))]; % Age
    TABLE1_VARS(2,groupi,:)=[mean(IDENTITY.DEMO(Sx,7)),std(IDENTITY.DEMO(Sx,7))]; % YrsEd
    for npi=1:5
        TABLE1_VARS(2+npi,groupi,:)=[nanmean(IDENTITY.NP(Sx,3+npi)),nanstd(IDENTITY.NP(Sx,3+npi))];
    end
end
CTL_Sx=logical( double(IDENTITY.DEMO(:,2)==1) .* double(IDENTITY.DEMO(:,3)==1)  );
Acute_Sx=logical( double(IDENTITY.DEMO(:,2)==1) .* double(IDENTITY.DEMO(:,3)==2)  );
Chronic_Sx=logical( double(IDENTITY.DEMO(:,2)==1) .* double(IDENTITY.DEMO(:,3)==3)  );

[~,T1_P,~,T1_STATS]=ttest2([IDENTITY.DEMO(CTL_Sx,[6,7]),IDENTITY.NP(CTL_Sx,4:8)],[IDENTITY.DEMO(Acute_Sx,[6,7]),IDENTITY.NP(Acute_Sx,4:8)])
[~,T1_P,~,T1_STATS]=ttest2([IDENTITY.DEMO(CTL_Sx,[6,7]),IDENTITY.NP(CTL_Sx,4:8)],[IDENTITY.DEMO(Chronic_Sx,[6,7]),IDENTITY.NP(Chronic_Sx,4:8)])

% TABLE 2
IDENTITY.TBI(Acute_Sx,4) % GCS
nanmedian(IDENTITY.TBI(Acute_Sx,6)) % LOCmins median
iqr(IDENTITY.TBI(Acute_Sx,6)) % LOCmins iqr
nansum(IDENTITY.TBI(Acute_Sx,7))
nanmedian(IDENTITY.TBI(Acute_Sx,8)) % Days median
iqr(IDENTITY.TBI(Acute_Sx,8)) % Days iqr

nanmedian(IDENTITY.TBI(Chronic_Sx,6)) % LOCmins median
iqr(IDENTITY.TBI(Chronic_Sx,6)) % LOCmins iqr
nansum(IDENTITY.TBI(Chronic_Sx,7))  % Data not here
nanmedian(IDENTITY.TBI(Chronic_Sx,9)) % Days median
iqr(IDENTITY.TBI(Chronic_Sx,9)) % Days iqr

%% ###########################################################################
clear INTERCOR_Rho INTERCOR_P;
% INDEX=[find(strcmp('BDI',IDENTITY_QUEX_HDR)),find(strcmp('NSItot',IDENTITY_QUEX_HDR)),find(strcmp('F_Tot_B4',IDENTITY_QUEX_HDR)),find(strcmp('F_Tot',IDENTITY_QUEX_HDR))]; 
COL={'bd','rd','md'};
SHIFT=[-.05,.05,0];
figure; hold on;
for groupi=1:2
    for timei=1:3
        Sx=logical( double(IDENTITY.DEMO(:,2)==timei) .* double(IDENTITY.DEMO(:,3)==groupi)  );
        for idxi=1:4
            subplot(2,2,idxi); hold on;
            ThisN=sum(~isnan(IDENTITY.QUEX(Sx,INDEX(idxi))));
            plot(timei+SHIFT(groupi),nanmean(IDENTITY.QUEX(Sx,INDEX(idxi))),COL{groupi});
            errorbar(timei+SHIFT(groupi),nanmean(IDENTITY.QUEX(Sx,INDEX(idxi))),nanstd(IDENTITY.QUEX(Sx,INDEX(idxi)))./sqrt(ThisN),'k.');
            set(gca,'xlim',[0 4],'xtick',[1:1:3]); 
            title(IDENTITY_QUEX_HDR{INDEX(idxi)});
            % ---
            Sx1=logical( double(IDENTITY.DEMO(:,2)==1) .* double(IDENTITY.DEMO(:,3)==groupi)  );
            Sx2=logical( double(IDENTITY.DEMO(:,2)==2) .* double(IDENTITY.DEMO(:,3)==groupi)  );
            Sx3=logical( double(IDENTITY.DEMO(:,2)==3) .* double(IDENTITY.DEMO(:,3)==groupi)  );
             plot([1+SHIFT(groupi) 2+SHIFT(groupi)],[nanmean(IDENTITY.QUEX(Sx1,INDEX(idxi))), nanmean(IDENTITY.QUEX(Sx2,INDEX(idxi)))],'k-');
             plot([2+SHIFT(groupi) 3+SHIFT(groupi)],[nanmean(IDENTITY.QUEX(Sx2,INDEX(idxi))), nanmean(IDENTITY.QUEX(Sx3,INDEX(idxi)))],'k-');
            % ---
        end
        if time==1
         [INTERCOR_Rho{groupi},INTERCOR_P{groupi}]=corr(IDENTITY.QUEX(Sx,[4,5,10]),'type','Spearman','rows','pairwise');
        end
    end
end
groupi=3; timei=1;
Sx=logical( double(IDENTITY.DEMO(:,2)==timei) .* double(IDENTITY.DEMO(:,3)==groupi)  );
for idxi=1:4
    subplot(2,2,idxi); hold on;
    ThisN=sum(~isnan(IDENTITY.QUEX(Sx,INDEX(idxi))));
    plot(timei+SHIFT(groupi),nanmean(IDENTITY.QUEX(Sx,INDEX(idxi))),COL{groupi});
    errorbar(timei+SHIFT(groupi),nanmean(IDENTITY.QUEX(Sx,INDEX(idxi))),nanstd(IDENTITY.QUEX(Sx,INDEX(idxi)))./sqrt(ThisN),'k.');
    set(gca,'xlim',[0 4],'xtick',[1:1:3]);
    title(IDENTITY_QUEX_HDR{INDEX(idxi)});
    
    
    [INTERCOR_Rho{3},INTERCOR_P{3}]=corr(IDENTITY.QUEX(Sx,[4,5,10]),'type','Spearman','rows','pairwise');
end


%%


