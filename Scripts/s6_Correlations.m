%%

COLS={'b','r','m'};
figure;
for si=1:3
    
    Sx=logical( double(IDENTITY.DEMO(:,2)==time) .* double(IDENTITY.DEMO(:,3)==si)  );
    
    % --------------------

    IV=squeeze(mean(MEGA_ERP(Sx,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)); 

    [rho,rho_p]=corr(IV,DV(Sx),'type','Spearman','rows','pairwise');
    [r,p]=corr(IV,DV(Sx),'type','Pearson','rows','pairwise');
    subplot(2,3,si); hold on
    scatter(IV,DV(Sx),COLS{si}); lsline
    set(gca,'xlim',[-10 20],'ylim',[20 120]);
 %    text(.1,.7,['df=',num2str(sum(logical(double(~isnan(IV)).*double(~isnan(DV(Sx))))) -2 ),' r=',num2str(r),' p=',num2str(p)],'sc');
    text(.1,.6,['df=',num2str(sum(logical(double(~isnan(IV)).*double(~isnan(DV(Sx))))) -2 ),' rho=',num2str(rho),' p=',num2str(rho_p)],'sc');
    clear rho rho_p r p IV
    title( BV_Chanlocs_60(ERPSITE(CONDI4Corr)).labels );
    
    % --------------------

    IV=squeeze(mean(MEGA_ERP(Sx,:,ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)); 
 
    [rho,rho_p]=corr(IV,DV(Sx),'type','Spearman','rows','pairwise');
    subplot(2,3,3+si); hold on
    rho_p(rho_p>=.05)=NaN; rho_p(rho_p<.05)=1; rho_p(isnan(rho_p))=0;
    topoplot(rho,BV_Chanlocs_60,'emarker2',{find(rho_p==1),'d','k',10,1});
    clear rho rhop_p IV;
    
    % --------------------
    
end

