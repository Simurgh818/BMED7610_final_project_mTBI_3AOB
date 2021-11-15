% mengz_JFC(r1, r2, r12, n) compares two correlations r1 and r2:
% r1: correlation between X and Y
% r2: correlation between X and Z
% r12: correlation between Y and Z
% n: number of observations used to compute correlations

%% ---------------
clear Sx CONDI4Corr rhoXY rhoXZ rhoYZ n menghyp mengp mengzscore

Sx=logical( double(IDENTITY.DEMO(:,2)==1) .* double(IDENTITY.DEMO(:,3)==3)  );
CONDI4Corr=2;
MING_CHRONIC.P3b=squeeze(mean(MEGA_ERP(Sx,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)); 
CONDI4Corr=3;
MING_CHRONIC.P3a=squeeze(mean(MEGA_ERP(Sx,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)); 

rhoXY=-.46;   % FrSBe & P3b
rhoXZ=.08;    % FrSBe & P3a
rhoYZ=corr(MING_CHRONIC.P3b,MING_CHRONIC.P3a,'type','Spearman','rows','pairwise');       % P3a & PBb
n=sum(Sx);
[menghyp,mengp,mengzscore] = mengz_JFC(rhoXY,rhoXZ,rhoYZ,n)

%% ---------------
clear Sx CONDI4Corr rhoXY rhoXZ rhoYZ n menghyp mengp mengzscore

Sx=logical( double(IDENTITY.DEMO(:,2)==1) .* double(IDENTITY.DEMO(:,3)==2)  );
CONDI4Corr=2;
MING_ACUTE.P3b=squeeze(mean(MEGA_ERP(Sx,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)); 
CONDI4Corr=3;
MING_ACUTE.P3a=squeeze(mean(MEGA_ERP(Sx,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)); 

rhoXY=-.11;   % FrSBe & P3b
rhoXZ=-.44;    % FrSBe & P3a
rhoYZ=corr(MING_ACUTE.P3b,MING_ACUTE.P3a,'type','Spearman','rows','pairwise');       % P3a & PBb
n=sum(Sx);
[menghyp,mengp,mengzscore] = mengz_JFC(rhoXY,rhoXZ,rhoYZ,n)

%% ---------------
clear Sx CONDI4Corr rhoXY rhoXZ rhoYZ n menghyp mengp mengzscore

Sx=logical( double(IDENTITY.DEMO(:,2)==2) .* double(IDENTITY.DEMO(:,3)==2)  );
CONDI4Corr=2;
MING_ACUTE_S2.P3b=squeeze(mean(MEGA_ERP(Sx,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)); 
CONDI4Corr=3;
MING_ACUTE_S2.P3a=squeeze(mean(MEGA_ERP(Sx,ERPSITE(CONDI4Corr),ERPWINS_tx2disp(CONDI4Corr,1):ERPWINS_tx2disp(CONDI4Corr,2),CONDI4Corr),3)); 

rhoXY=-.11;   % FrSBe & P3b
rhoXZ=-.49;    % FrSBe & P3a
rhoYZ=corr(MING_ACUTE_S2.P3b,MING_ACUTE_S2.P3a,'type','Spearman','rows','pairwise');       % P3a & PBb
n=sum(Sx);
[menghyp,mengp,mengzscore] = mengz_JFC(rhoXY,rhoXZ,rhoYZ,n)

