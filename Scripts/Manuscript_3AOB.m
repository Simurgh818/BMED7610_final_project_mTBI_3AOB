%% Step 3 3AOB
clear all; clc
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc');
datadir='C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\';
homedir='C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\Scripts';
cd(homedir);

load('C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\Scripts\BV_Chanlocs_60.mat');
tx2disp=-500:2:1000;

% Load Data
s1_Load_Data;       

% Kill Data
s2_Kill_Data;

%% Demographics

% s3_Demographics;
 
% sx_Predict_Attrition;


%% Example ERPs

StdSite=find(strcmpi('FCz',{BV_Chanlocs_60.labels}));
StdT1=300; StdT2=450; 

TargSite=find(strcmpi('Pz',{BV_Chanlocs_60.labels}));
TargT1=400; TargT2=600; 

NovSite=find(strcmpi('FCz',{BV_Chanlocs_60.labels}));
NovT1=300; NovT2=450; 

ERPSITE=[StdSite,TargSite,NovSite];
ERPWINS=[StdT1,StdT2;TargT1,TargT2;NovT1,NovT2];
ERPWINS_tx2disp=[[find(tx2disp==StdT1),find(tx2disp==StdT2)];...
                 [find(tx2disp==TargT1),find(tx2disp==TargT2)];...
                 [find(tx2disp==NovT1),find(tx2disp==NovT2)] ];


s4_Example_ERPs

%% ERPs by Group

time=1;

s5_ERPs_by_Group


%% For SPSS

s6_FOR_SPSS_s2

figure; boxplot(FORSPSS(:,[10,11,16,17,22,23])); % Raw, Scaled
skewness(FORSPSS(:,[10,11,16,17,22,23]))  % not skewed

% Calculate reliability for controls
CTL_REL=FORSPSS(FORSPSS(:,3)==1,:);

[REL.rho.F12,REL.p.F12]=corr(CTL_REL(:,11),CTL_REL(:,17),'type','Spearman','rows','pairwise'); % F_Tot 1 & 2
[REL.rho.F13,REL.p.F13]=corr(CTL_REL(:,11),CTL_REL(:,23),'type','Spearman','rows','pairwise'); % F_Tot 1 & 3
[REL.rho.F23,REL.p.F23]=corr(CTL_REL(:,17),CTL_REL(:,23),'type','Spearman','rows','pairwise'); % F_Tot 2 & 3

[REL.rho.P3b12,REL.p.P3b12]=corr(CTL_REL(:,11+1),CTL_REL(:,17+1),'type','Spearman','rows','pairwise'); % P3b 1 & 2
[REL.rho.P3b13,REL.p.P3b13]=corr(CTL_REL(:,11+1),CTL_REL(:,23+1),'type','Spearman','rows','pairwise'); % P3b 1 & 3
[REL.rho.P3b23,REL.p.P3b23]=corr(CTL_REL(:,17+1),CTL_REL(:,23+1),'type','Spearman','rows','pairwise'); % P3b 2 & 3

[REL.rho.P3a12,REL.p.P3a12]=corr(CTL_REL(:,11+2),CTL_REL(:,17+2),'type','Spearman','rows','pairwise'); % P3a 1 & 2
[REL.rho.P3a13,REL.p.P3a13]=corr(CTL_REL(:,11+2),CTL_REL(:,23+2),'type','Spearman','rows','pairwise'); % P3a 1 & 3
[REL.rho.P3a23,REL.p.P3a23]=corr(CTL_REL(:,17+2),CTL_REL(:,23+2),'type','Spearman','rows','pairwise'); % P3a 2 & 3

%% Correlations

DV=IDENTITY.QUEX(:,find(strcmp('F_Tot',IDENTITY_QUEX_HDR)));  

time=1; 
CONDI4Corr=3;  % Std, Targ, Nov

s6_Correlations

%% Predictions

quexidx=find(strcmp('F_Tot',IDENTITY_QUEX_HDR));
CONDI4Corr=2;  % Std, Targ, Nov

s6_Correlations_S1EEG_With_FrSBediffs

%% -------------- Between-Group rho-to-z

% Just type 'em in here from the plots (remember number on plots is df, not N):

r1=-.11
n1=38
r2=-.46
n2=23

clc;

t_r1 = 0.5*log((1+r1)/(1-r1));
t_r2 = 0.5*log((1+r2)/(1-r2));
z = (t_r1-t_r2)/sqrt(1/(n1-3)+1/(n2-3))
p = (1-normcdf(abs(z),0,1))*2


%% -------------- Within-Group rho-to-z

s7_Mengs_z


