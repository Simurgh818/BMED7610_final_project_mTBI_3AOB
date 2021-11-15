%%

TITLES={'Std','Targ','Nov'};
COL={'b','r','m'};
YLIM=[-6 8];   


figure;
for ai=1:3
    subplot(3,1,ai); hold on;
     for gi=1:3
        V = logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==gi) );
        plot(tx2disp,squeeze(nanmean( MEGA_ERP(V,ERPSITE(ai),:,ai) ,1)),COL{gi},'Linewidth',2);  
        clear V;
    end
    for gi=1:3
        V = logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==gi) );
        shadedErrorBar(tx2disp,squeeze(nanmean( MEGA_ERP(V,ERPSITE(ai),:,ai) ,1)),squeeze(nanstd( MEGA_ERP(V,ERPSITE(ai),:,ai) ,1))./sqrt(sum(V)),COL{gi});
        clear V;
    end
    for gi=1:3
        V = logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==gi) );
        plot(tx2disp,squeeze(nanmean( MEGA_ERP(V,ERPSITE(ai),:,ai) ,1)),COL{gi},'Linewidth',2);   
        set(gca,'ylim',YLIM);
        clear V;
    end
    
    plot([ERPWINS(ai,1) ERPWINS(ai,1)],[5 7],'k:');   plot([ERPWINS(ai,2) ERPWINS(ai,2)],[5 7],'k:'); 
    
    title(TITLES{ai});
    
    V_ctl = logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==1) );
    V_acute = logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==2) );
    V_chronic = logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==3) );
    ERPs_ctl=squeeze(MEGA_ERP(V_ctl,ERPSITE(ai),:,ai));
    ERPs_acute=squeeze(MEGA_ERP(V_acute,ERPSITE(ai),:,ai));
    ERPs_chronic=squeeze(MEGA_ERP(V_chronic,ERPSITE(ai),:,ai));
    [H,P,CI,STATS]=ttest2(ERPs_ctl,ERPs_acute); P(P>.05)=NaN;  P(P<=.05)=1; 
    plot(tx2disp,-4.*P,'r');  clear H P CI STATS TEMP*;
    [H,P,CI,STATS]=ttest2(ERPs_ctl,ERPs_chronic); P(P>.05)=NaN;  P(P<=.05)=1; 
    plot(tx2disp,-4.5.*P,'m');  clear H P CI STATS TEMP*;    
    
    clear V_ctl V_acute V_chronic ERPs_ctl ERPs_acute ERPs_chronic
     
end
legend({'CTL','Acute','Chronic'},'Location','NorthWest')


%%

figure;
for ai=1:3  % Condi
    subplot(3,1,ai); hold on;
    for gi=1:3 % Group
        for time=1:3
            V = logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==gi) );
            plot(time,squeeze(nanmean(mean( MEGA_ERP(V,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3),1)),[COL{gi},'d']);
            errorbar(time,squeeze(nanmean(mean( MEGA_ERP(V,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3),1)),...
                          squeeze(nanstd(mean( MEGA_ERP(V,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3),1)) ./sqrt(sum(V)),'k.');     
            BIG_OL_N(gi,time)=sum(V);       
                        
            clear V;
        end
        
        % ---
        V1 = logical( double(IDENTITY.DEMO(:,2)==1) .* double(IDENTITY.DEMO(:,3)==gi) );
        V2 = logical( double(IDENTITY.DEMO(:,2)==2) .* double(IDENTITY.DEMO(:,3)==gi) );
        V3 = logical( double(IDENTITY.DEMO(:,2)==3) .* double(IDENTITY.DEMO(:,3)==gi) );
        plot([1 2],[squeeze(nanmean(mean( MEGA_ERP(V1,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3),1)) ,...
            squeeze(nanmean(mean( MEGA_ERP(V2,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3),1))],'k-');
        plot([2 3],[squeeze(nanmean(mean( MEGA_ERP(V2,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3),1)) ,...
            squeeze(nanmean(mean( MEGA_ERP(V3,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3),1))],'k-');
        clear V*;
        % ---
        
    end

    set(gca,'xlim',[0 4],'xtick',[1:1:3])
    title(TITLES{ai});
end

%%

STATS{1}=[]; STATS{2}=[];  % Targ, Nov
for si=1:2  % CTL, Acute
    Sx=logical( double(IDENTITY.DEMO(:,3)==si)  );
    Sx_idxs=unique(IDENTITY.DEMO(Sx,1));
    for sxi=1:length(Sx_idxs)
        thisguy=Sx_idxs(sxi);
        
        FIRST=[]; SECOND=[]; THIRD=[];
        
        FIRST=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==1) ));
        SECOND=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==2) ));
        THIRD=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==3) ));
        
        ai=2;
        ERP1=NaN; ERP2=NaN; ERP3=NaN;
        if ~isempty(FIRST), ERP1=mean( MEGA_ERP(FIRST,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3); end
        if ~isempty(SECOND), ERP2=mean( MEGA_ERP(SECOND,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3); end
        if ~isempty(THIRD), ERP3=mean( MEGA_ERP(THIRD,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3); end
        STATS{ai-1}=[STATS{ai-1};thisguy,si,ERP1,ERP2,ERP3];

        ai=3;
        ERP1=NaN; ERP2=NaN; ERP3=NaN;
        if ~isempty(FIRST), ERP1=mean( MEGA_ERP(FIRST,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3); end
        if ~isempty(SECOND), ERP2=mean( MEGA_ERP(SECOND,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3); end
        if ~isempty(THIRD), ERP3=mean( MEGA_ERP(THIRD,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3); end
        STATS{ai-1}=[STATS{ai-1};thisguy,si,ERP1,ERP2,ERP3];
        
        clear thisguy;
    end
    clear Sx Sx_idxs;
end

si=3; clear Sx Sx_idxs; % chronic
Sx=logical( double(IDENTITY.DEMO(:,3)==si)  );
Sx_idxs=unique(IDENTITY.DEMO(Sx,1));
for sxi=1:length(Sx_idxs)
    
    thisguy=Sx_idxs(sxi);
    FIRST=[];
    FIRST=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==1) ));
    
    ai=2;
    ERP1=NaN;
    if ~isempty(FIRST), ERP1=mean( MEGA_ERP(FIRST,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3); end
    STATS{ai-1}=[STATS{ai-1};thisguy,si,ERP1,NaN,NaN];
    
    ai=3;
    ERP1=NaN;
    if ~isempty(FIRST), ERP1=mean( MEGA_ERP(FIRST,ERPSITE(ai),ERPWINS_tx2disp(ai,1):ERPWINS_tx2disp(ai,2),ai) ,3); end
    STATS{ai-1}=[STATS{ai-1};thisguy,si,ERP1,NaN,NaN];
    
    clear thisguy;
end

