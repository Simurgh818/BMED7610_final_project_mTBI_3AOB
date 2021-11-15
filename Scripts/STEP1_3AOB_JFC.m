%% 3AOB JFC
clear
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\');
addpath(genpath('C:\Users\sinad\eeglab2021.1')); 
rmpath('C:\Users\sinad\eeglab2021.1\functions');
rmpath('C:\Users\sinad\eeglab2021.1\functions');
datadir='C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\'; % Data are here
saveloc='C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\';
load('C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\Scripts\BV_Chanlocs_60.mat');
cd(saveloc);

sx_dirs=dir([datadir,'*.mat']); 
for sxi=1:length(sx_dirs)
    for ses=1:3
        sessdir=[datadir,sx_dirs(sxi).name,'C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\ds003522\sub-001\ses-01\eeg\'];
        sx_sess{sxi}{1,ses}=dir([sessdir,'*_',num2str(ses),'_ODDBALL.vhdr']);
        sx_sess{sxi}{2,ses}=sessdir;
        sx_sess{sxi}{3,ses}=sx_dirs(sxi).name;
    end
end
LOG=[];
for sxi=1:size(sx_sess,2)
    for sess=1:3
        if ~isempty( sx_sess{sxi}{1,sess} )
            subno=str2num(sx_sess{sxi}{1,sess}.name(1:4));
            URSI=sx_sess{sxi}{3,sess};
            LOG(subno-3000,sess+1)=subno;
            LOG(subno-3000,1)=str2num(URSI(end-4:end));
        end
    end
end

for sxi=1:size(sx_sess,2)
    for sess=1:3
        if ~isempty( sx_sess{sxi}{1,sess} )
            
            subno=str2num(sx_sess{sxi}{1,sess}.name(1:4));
            thisdir=sx_sess{sxi}{2,sess};
            URSI=sx_sess{sxi}{3,sess};
            LOG2(sxi,sess)=subno;
            
            if ~exist([saveloc,num2str(subno),'_',num2str(sess),'_3AOB.mat']);
                
                % Data are 65 chans: 1=63 is EEG, 64 is VEOG, 65 is EKG  Ref'd to CPz - - will want to retrieve that during re-referencing
                EEG = pop_loadbv(thisdir,[num2str(subno),'_',num2str(sess),'_ODDBALL.vhdr']);  clc; disp(['Loading ',num2str(subno),'  s',num2str(sess)]);
                % Run PATCH for sx<3003 s<2   AND    for bad templates
                PATCH
                % Get Locs
                locpath=('Y:\Programs\eeglab12_0_2_1b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp');
                EEG = pop_chanedit(EEG,    'lookup', locpath);
                EEG = eeg_checkset( EEG );
                % Get event types
                for ai=2:length(EEG.event); clear temp; temp=EEG.event(ai).type;
                    if isempty(strmatch('boundary',temp)); TYPES(ai)=str2num(temp(2:end)) ; clear temp; end
                end
                UNIQUE_TYPES=unique(TYPES);
                for ai=1:length(UNIQUE_TYPES); UNIQUE_TYPES_COUNT(ai)=sum(TYPES==UNIQUE_TYPES(ai)); end
                clc; TRIGGERS=[UNIQUE_TYPES;UNIQUE_TYPES_COUNT] % Trigger type, Frequency
                
                % Epoch
                All_STIM={'S201','S200','S202'};   % Std, Target, Novel
                EEG = pop_epoch( EEG, All_STIM, [-2 2], 'newname', 'Epochs', 'epochinfo', 'yes');
                EEG = eeg_checkset( EEG );
                % Remove VEOG and EKG
                EEG.EKG=squeeze(EEG.data(65,:,:));
                EEG.VEOG=squeeze(EEG.data(64,:,:));
                EEG.data=EEG.data(1:63,:,:);
                EEG.nbchan=63;
                EEG.chanlocs(65)=[]; EEG.chanlocs(64)=[];
                % Fix BV-specific issue - - - only needed for APPLE
                for ai=1:size(EEG.urevent,2), EEG.urevent(ai).bvtime=EEG.urevent(ai).bvmknum; end
                for ai=1:size(EEG.event,2), EEG.event(ai).bvtime=EEG.event(ai).bvmknum; end
                for ai=1:size(EEG.epoch,2), EEG.epoch(ai).eventbvtime=EEG.epoch(ai).eventbvmknum; end
                % Add CPz
                EEG = pop_chanedit(EEG,'append',63,'changefield',{64 'labels' 'CPz'});
                EEG = pop_chanedit(EEG,'lookup', locpath);
                % Re-Ref to Average Ref and recover CPz
                EEG = pop_reref(EEG,[],'refloc',struct('labels',{'CPz'},'type',{''},'theta',{180},'radius',{0.12662},'X',{-32.9279},'Y',{-4.0325e-15},'Z',{78.363},...
                    'sph_theta',{-180},'sph_phi',{67.208},'sph_radius',{85},'urchan',{64},'ref',{''}),'keepref','on');
                % Remove everything else NOW that CPz has been reconstructed from the total
                EEG.MASTOIDS = squeeze(mean(EEG.data([10,21],:,:),1));
                EEG.data = EEG.data([1:4,6:9,11:20,22:26,28:64],:,:);
                EEG.nbchan=60;
                EEG.chanlocs(27)=[];  EEG.chanlocs(21)=[];   EEG.chanlocs(10)=[];   EEG.chanlocs(5)=[];  % Have to be in this order!
                % Should probably re-ref to average again now that the contaminated channels are gone
                EEG = pop_reref(EEG,[]);
                % Remove mean
                EEG = pop_rmbase(EEG,[],[]);
                
                
                % ----------------------
                % Setup APPLE to interp chans, reject epochs, & ID bad ICs.  Output will be Avg ref'd and ICA'd.
                eeg_chans=1:60;
                Do_ICA=1;
                ref_chan=36;  % Re-Ref to FCz   [WEIRD STEP, BUT THIS IS FOR FASTER, which is a part of APPLE]
                EEG = pop_reref(EEG,ref_chan,'keepref','on');
                
                % Run APPLE (will re-ref data to avg ref)
                [EEG,EEG.bad_chans,EEG.bad_epochs,EEG.bad_ICAs]=APPLE_3AOB(EEG,eeg_chans,ref_chan,Do_ICA,subno,EEG.VEOG,sess,BV_Chanlocs_60);
                
                % Save
                save([num2str(subno),'_',num2str(sess),'_3AOB.mat'],'EEG');
                % ----------------------
                
                %% Remove the (presumptive) bad ICAs:
                bad_ICAs_To_Remove=EEG.bad_ICAs{2};
                if bad_ICAs_To_Remove==0, bad_ICAs_To_Remove=1; end
                EEG = pop_subcomp( EEG, bad_ICAs_To_Remove, 0);
                
                
                % Get the good info out of the epochs
                for ai=1:size(EEG.epoch,2)
                    % Initialize
                    EEG.epoch(ai).CUE=NaN;
                    for bi=1:size(EEG.epoch(ai).eventlatency,2)
                        % Get STIMTYPE
                        if EEG.epoch(ai).eventlatency{bi}==0 && isempty(strmatch(EEG.epoch(ai).eventtype{bi},'N999')); % If this bi is the event
                            % Get StimType
                            FullName=EEG.epoch(ai).eventtype{bi};
                            EEG.epoch(ai).CUE=str2num(FullName(2:end)) ;
                            clear FullName
                            VECTOR(ai,1)=EEG.epoch(ai).CUE;
                        end
                    end
                end
                
                % Let's just do this for display
                dims=size(EEG.data);
                EEG.data=eegfilt(EEG.data,500,[],20);
                EEG.data=reshape(EEG.data,dims(1),dims(2),dims(3));
                
                % Set Params
                tx=-2000:2:1998;
                b1=find(tx==-200);  b2=find(tx==0);
                t1=find(tx==-500);  t2=find(tx==1000);
                toporange1=find(tx==250);  toporange2=find(tx==600);  toporangetot=250:2:600;
                tx2disp=-500:2:1000;
                MAPLIMS=[-8 8];
                
                % Basecor your ERPs here so they are pretty.
                BASE=squeeze(  mean(EEG.data(:,b1:b2,:),2)  );
                for ai=1:dims(1)
                    EEG.data(ai,:,:)=squeeze(EEG.data(ai,:,:))-repmat( BASE(ai,:),dims(2),1 );
                end
                
                
                % Get max of P2 across all condis
                site=11;  % Pz
                ERP4topo=mean(EEG.data(site,toporange1:toporange2,VECTOR(:,1)==200),3);
                topomax_P3b=toporangetot(find(ERP4topo==max(ERP4topo)));
                topotoplot_P3b=find(tx==topomax_P3b);
                site=36;  % FCz
                ERP4topo=mean(EEG.data(site,toporange1:toporange2,VECTOR(:,1)==202),3);
                topomax_P3a=toporangetot(find(ERP4topo==max(ERP4topo)));
                topotoplot_P3a=find(tx==topomax_P3a);
                % --------------
                figure;
                site=11;  % Pz
                subplot(3,4,1:4); hold on
                plot(tx2disp,mean(EEG.data(site,t1:t2,VECTOR(:,1)==201),3),'k');
                plot(tx2disp,mean(EEG.data(site,t1:t2,VECTOR(:,1)==200),3),'r');
                plot(tx2disp,mean(EEG.data(site,t1:t2,VECTOR(:,1)==202),3),'b');
                plot([topomax_P3b topomax_P3b],[-2 2],'m','linewidth',2); % indicate the max with a magenta line
                title(['Pz   Subno: ',num2str(subno),'  Sess:',num2str(sess)]);
                legend({'Std','Target','Novel'},'Location','NorthWest');
                % --------------
                site=36;  % FCz
                subplot(3,4,5:8); hold on
                plot(tx2disp,mean(EEG.data(site,t1:t2,VECTOR(:,1)==201),3),'k');
                plot(tx2disp,mean(EEG.data(site,t1:t2,VECTOR(:,1)==200),3),'r');
                plot(tx2disp,mean(EEG.data(site,t1:t2,VECTOR(:,1)==202),3),'b');
                plot([topomax_P3a topomax_P3a],[-2 2],'m','linewidth',2); % indicate the max with a magenta line
                title(['FCz   Subno: ',num2str(subno),'  Sess:',num2str(sess)]);
                % --------------
                subplot(3,4,9); hold on
                topoplot( mean(EEG.data(:,topotoplot_P3b,VECTOR(:,1)==201),3) , BV_Chanlocs_60,'maplimits',MAPLIMS); title('Std @ P3b')
                subplot(3,4,10); hold on
                topoplot( mean(EEG.data(:,topotoplot_P3b,VECTOR(:,1)==200),3) , BV_Chanlocs_60,'maplimits',MAPLIMS); title('Targ')
                subplot(3,4,11); hold on
                topoplot( mean(EEG.data(:,topotoplot_P3a,VECTOR(:,1)==202),3) , BV_Chanlocs_60,'maplimits',MAPLIMS); title('Novel')
                
                saveas(gcf, [num2str(subno),'_',num2str(sess),'_3AOB_ERPs.png'],'png');
                close all;
                
                clear EEG VECTOR BASE PROBE TRIGGERS TYPES UNIQUE* did* topo* ERP* URSI dims eeg_chans;
            end
        end
    end
end

%%
 
