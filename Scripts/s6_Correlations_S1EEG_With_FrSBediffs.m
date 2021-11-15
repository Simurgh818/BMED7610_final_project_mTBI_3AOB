%%
COLS={'c','m'};
ROWS={'1-2','1-3','2-3'};
figure;
for si=1:2
    Sx=logical( double(IDENTITY.DEMO(:,3)==si)  );
    Sx_idxs=unique(IDENTITY.DEMO(Sx,1));
    
    clear IV* DV* *12 *23 *13; 
    for sxi=1:length(Sx_idxs)
        thisguy=Sx_idxs(sxi);

        FIRST=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==1) ));
        SECOND=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==2) ));
        THIRD=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==3) ));
        
        if ~isempty(FIRST) && ~isempty(SECOND)
            IVs12(sxi,:)= squeeze(mean(MEGA_ERP(FIRST,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3))  ;
            DVs12(sxi,:)=IDENTITY.QUEX(SECOND,quexidx)-IDENTITY.QUEX(FIRST,quexidx);
            age12(sxi,:)=IDENTITY.QUEX(FIRST,6);
            TOPF12(sxi,:)=IDENTITY.NP(FIRST,4);
            sex12(sxi,:)=IDENTITY.DEMO(FIRST,5);
            % NP vars that also predicted dropout
            Span12(sxi,:)=IDENTITY.NP(FIRST,6);
            Coding12(sxi,:)=IDENTITY.NP(FIRST,5);
        else
            IVs12(sxi,:)=NaN;
            DVs12(sxi,:)=NaN;
            age12(sxi,:)=NaN;
            TOPF12(sxi,:)=NaN;
            sex12(sxi,:)=NaN;
        end
        
       if ~isempty(SECOND) && ~isempty(THIRD)
            IVs23(sxi,:)= squeeze(mean(MEGA_ERP(SECOND,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)) ;
            DVs23(sxi,:)=IDENTITY.QUEX(THIRD,quexidx)-IDENTITY.QUEX(SECOND,quexidx);
        else
            IVs23(sxi,:)=NaN;
            DVs23(sxi,:)=NaN;
        end
        
        if ~isempty(FIRST) && ~isempty(THIRD)
            IVs13(sxi,:)= squeeze(mean(MEGA_ERP(FIRST,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3))  ;
            DVs13(sxi,:)=IDENTITY.QUEX(THIRD,quexidx)-IDENTITY.QUEX(FIRST,quexidx);
        else
            IVs13(sxi,:)=NaN;
            DVs13(sxi,:)=NaN;
        end

    end

    % --------------------

    [rho,rho_p]=corr(IVs12,DVs12,'type','Spearman','rows','pairwise');
    [r,p]=corr(IVs12,DVs12,'type','Pearson','rows','pairwise');
    subplot(3,2,si); hold on
    scatter(IVs12,DVs12,COLS{si}); lsline
    set(gca,'xlim',[-10 20],'ylim',[-40 40]);
    text(.1,.7,['df=',num2str(sum(logical(double(~isnan(IVs12)).*double(~isnan(DVs12)))) -2 ),' r=',num2str(r),' p=',num2str(p)],'sc');
    text(.1,.6,['df=',num2str(sum(logical(double(~isnan(IVs12)).*double(~isnan(DVs12)))) -2 ),' rho=',num2str(rho),' p=',num2str(rho_p)],'sc');
    clear rho rho_p r p IV
    title( [BV_Chanlocs_60(ERPSITE(CONDI4Corr)).labels,'   S1EEG    ',ROWS{1}] );
    
    % --------------------

    [rho,rho_p]=corr(IVs13,DVs13,'type','Spearman','rows','pairwise');
    [r,p]=corr(IVs13,DVs13,'type','Pearson','rows','pairwise');
    subplot(3,2,si+2); hold on
    scatter(IVs13,DVs13,COLS{si}); lsline
    set(gca,'xlim',[-10 20],'ylim',[-40 40]);
    text(.1,.7,['df=',num2str(sum(logical(double(~isnan(IVs13)).*double(~isnan(DVs13))))  -2 ),' r=',num2str(r),' p=',num2str(p)],'sc');
    text(.1,.6,['df=',num2str(sum(logical(double(~isnan(IVs13)).*double(~isnan(DVs13))))  -2 ),' rho=',num2str(rho),' p=',num2str(rho_p)],'sc');
    clear rho rho_p r p IV
    title( [BV_Chanlocs_60(ERPSITE(CONDI4Corr)).labels,'   S1EEG    ',ROWS{2}] );
    
    % --------------------

    [rho,rho_p]=corr(IVs23,DVs23,'type','Spearman','rows','pairwise');
    [r,p]=corr(IVs23,DVs23,'type','Pearson','rows','pairwise');
    subplot(3,2,si+4); hold on
    scatter(IVs23,DVs23,COLS{si}); lsline
    set(gca,'xlim',[-10 20],'ylim',[-40 40]);
    text(.1,.7,['df=',num2str(sum(logical(double(~isnan(IVs23)).*double(~isnan(DVs23))))  -2 ),' r=',num2str(r),' p=',num2str(p)],'sc');
    text(.1,.6,['df=',num2str(sum(logical(double(~isnan(IVs23)).*double(~isnan(DVs23))))  -2 ),' rho=',num2str(rho),' p=',num2str(rho_p)],'sc');
    clear rho rho_p r p IV
    title( [BV_Chanlocs_60(ERPSITE(CONDI4Corr)).labels,'   S1EEG    ',ROWS{3}] );
     
    % --------------------

end

 
% Check demographic (S1) vars in the sub-acute group on FrSBe change

[r,p]=corr(DVs12,age12,'type','Spearman','rows','pairwise')

[r,p]=corr(DVs12,TOPF12,'type','Spearman','rows','pairwise')

[H,P,CI,STATS]=ttest2(DVs12(sex12==1),DVs12(sex12==0))

[r,p]=corr(DVs12,Span12,'type','Spearman','rows','pairwise')
[r,p]=corr(DVs12,Coding12,'type','Spearman','rows','pairwise')

   
   

[r,p]=corr(DVs13,age12,'type','Spearman','rows','pairwise')
[r,p]=corr(DVs13,TOPF12,'type','Spearman','rows','pairwise')
[H,P,CI,STATS]=ttest2(DVs13(sex12==1),DVs13(sex12==0))
[r,p]=corr(DVs13,Span12,'type','Spearman','rows','pairwise')
[r,p]=corr(DVs13,Coding12,'type','Spearman','rows','pairwise')
