%%

acuteTBI=IDENTITY.DEMO(:,3)==2;
Sx_idxs=unique(IDENTITY.DEMO(acuteTBI,1));

AttritionPredictors=NaN(length(Sx_idxs),14);
for sxi=1:length(Sx_idxs)
    thisguy=Sx_idxs(sxi);
    
    FIRST=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==1) ));
    SECOND=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==2) ));
    THIRD=find(logical( double(IDENTITY.DEMO(:,1)==thisguy) .* double(IDENTITY.DEMO(:,2)==3) ));
    
    if ~isempty(FIRST),  Attrition(sxi,1)=1;
        AttritionPredictors_HDR={'age';'sex';'TOPF';'Coding';'Span';'HVLT13';'HVLTDelay';'GCS';'LOCtime';'LOM';'Days';'BDI';'NSI';'FrSBe'};
        AttritionPredictors(sxi,:)=[...
            IDENTITY.DEMO(FIRST,6:7),...
            IDENTITY.NP(FIRST,4:8),...
            IDENTITY.TBI(FIRST,[4,6:8]),...
            IDENTITY.QUEX(FIRST,[4,5,10]),...
            ];
    else  Attrition(sxi,1)=0; end
    if ~isempty(SECOND),  Attrition(sxi,2)=1; else  Attrition(sxi,2)=0; end
    if ~isempty(THIRD),  Attrition(sxi,3)=1; else  Attrition(sxi,3)=0; end
    
    clear thisguy FIRST SECOND THIRD;
end

% Clear NaNs
AttritionPredictors=AttritionPredictors(~isnan(AttritionPredictors(:,1)),:);
Attrition=Attrition(~isnan(AttritionPredictors(:,1)),:);


for atti=1:size(AttritionPredictors,2)
    for sessi=1:3  % 1 doesn't really make sense, but added here to keep columns nice
        
        A=AttritionPredictors(:,atti);
        B=Attrition(:,sessi);
        
        A2=A(~isnan(A));
        B2=B(~isnan(A));
        
        % B_acc is [constant, v1, v2,v1*v2] [validated by SPSS]
        % STATS_acc.p is the p value for each
        [B_acc,DEV_acc,STATS_acc] = glmfit(zscore(A2),B2, 'binomial','link','logit');
        Predictors_Logistic{sessi}(atti,1)=B_acc(2);
        Predictors_Logistic{sessi}(atti,2)=STATS_acc.p(2);
        
        [~,P,~,STATS]=ttest2(A2(B2==1),A2(B2==0));
        Predictors_t{sessi}(atti,1)=STATS.tstat;
        Predictors_t{sessi}(atti,2)=P;
        
        [P,~,U]=ranksum(A2(B2==1),A2(B2==0));
        Predictors_u{sessi}(atti,1)=U.zval;
        Predictors_u{sessi}(atti,2)=P;
    
         clear B_acc DEV_acc STATS_acc A B A2 B2 STATS P U;
    end
end

% % % S2: 
% Logistic - 5
% t-test - 5
% U-test - 5

% % % S3: 
% Logistic - 3,5 [almost 4]
% t-test - 3,5  [almost 4]
% U-test - 3,4,5

% 3=TOPF
% 4=Coding
% 5=Span

% Sess 2 dropouts:
nanmean(AttritionPredictors(Attrition(:,2)==0,3:5))
% Sess 2 stays:
nanmean(AttritionPredictors(Attrition(:,2)==1,3:5))

% Sess 3 dropouts:
nanmean(AttritionPredictors(Attrition(:,3)==0,3:5))
% Sess 3 stays:
nanmean(AttritionPredictors(Attrition(:,3)==1,3:5))

% So Lower Span predicts S2 dropout and Lower Span, Coding, and TOPF predict S3 dropout
