%% Step 2 Oddball

clear all; clc
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\');
savedir='C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\';

load('C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\Scripts\BV_Chanlocs_60.mat');

% ########## For Cavanagh data
datadir='C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\';
[D_DAT,D_HDR,D_ALL]=xlsread('C:\Users\sinad\OneDrive - Georgia Institute of Technology\BMED7610QuantitativeNeuro-CompNeuro\Final Project\mTBI 3AOB PreProc\BMED7610_final_project_mTBI_3AOB\Scripts\QUALITY_CHECK.xlsx','ODDBALL_ICAs');
FILEENDER='_3AOB.mat';

% % % ########## For Quinn data
% % datadir='Z:\EXPERIMENTS\mTBICoBRE\EEG\QUINN 3AOB Preproc\';
% % [D_DAT,D_HDR,D_ALL]=xlsread('Z:\EXPERIMENTS\mTBICoBRE\ANALYSIS\QUINN_QUALITY_CHECK.xlsx','ODDBALL_ICAs');
% % FILEENDER='_QUINN_3AOB.mat';

cd(datadir);

% ############# Set Params
srate=500;
tx=-2000:1000/srate:1998;
B1=find(tx==-300);  B2=find(tx==-200);
T1=find(tx==-500);  T2=find(tx==1000);
tx2disp=-500:2:1000;
% #############


for si=1:length(D_DAT)
    for sess=1:size(D_DAT,2)-1 % should be '2' for Quinn data, '3' for Cavanagh data
        
        subno=D_DAT(si,1);
        skip=0;
        
        INFO=D_ALL{si+1,sess+1};  % +1's b/c of subno column and header row
        disp(['TRYOUT ',num2str(subno),' S',num2str(sess)]);
        
        if isnumeric(INFO), bad_ICAs_To_Remove=INFO; end
        if isnan(INFO), skip=1; end     % not done yet
        if strmatch('BAD',INFO), skip=1; end  % Bad data
        if ~isnumeric(INFO), bad_ICAs_To_Remove=str2num(INFO); end
        
        % Don't repeat if already done
        if exist([savedir,num2str(subno),'_',num2str(sess),'_3AOB_TFandERPs_L.mat'])==2, skip=1; end
        
        if skip==0
            
            load([num2str(subno),'_',num2str(sess),FILEENDER]);  disp(['DOING: ',num2str(subno),'_',num2str(sess),'_3AOB.mat']);
            
            % Remove the bad ICAs:
            disp(['BAD ICAS:  ', num2str(bad_ICAs_To_Remove)]);
            EEG = pop_subcomp( EEG, bad_ICAs_To_Remove, 0);

           % Get the good info out of the epochs
            for ai=1:size(EEG.epoch,2)
                % Initialize
                EEG.epoch(ai).EEG=NaN;
                for bi=1:size(EEG.epoch(ai).eventlatency,2)
                    % Get STIMTYPE
                    if EEG.epoch(ai).eventlatency{bi}==0 && isempty(strmatch(EEG.epoch(ai).eventtype{bi},'N999')); % If this bi is the event
                        % Get StimType
                        FullName=EEG.epoch(ai).eventtype{bi};
                        EEG.epoch(ai).EEG=str2num(FullName(2:end)) ;
                        
                        clear FullName
                        VECTOR(ai,1)=EEG.epoch(ai).EEG;  All_STIM={'S201','S200','S202'};   % Std, Target, Novel
                    end
                end
            end
            
            % Only as many STD as NOV
            N_n=sum(VECTOR(:,1)==202);
            temp_idxs=find(VECTOR(:,1)==201);
            temp_idxs=shuffle(temp_idxs);
            VECTOR(temp_idxs(N_n+1:end),1)=999;  clear temp_idxs;
            % Save trial counts
            TRL_ct(1)=sum(VECTOR(:,1)==201);
            TRL_ct(2)=sum(VECTOR(:,1)==200);
            TRL_ct(3)=sum(VECTOR(:,1)==202);
      
            
            %%
            % $$$$$$$$$$$$$$$$$$$$$$$           $$$$$$$$$$$$$$$$$$$$$$$
            % $$$$$$$$$$$$$$$$$$$$$$$ Time-Freq
            % $$$$$$$$$$$$$$$$$$$$$$$           $$$$$$$$$$$$$$$$$$$$$$$
            
            % Setup Wavelet Params
            num_freqs=50;
            frex=logspace(.01,1.7,num_freqs);
            s=logspace(log10(3),log10(10),num_freqs)./(2*pi*frex);
            t=-2:1/EEG.srate:2;
            
            % Definte Convolution Parameters
            dims = size(EEG.data);
            n_wavelet = length(t);
            n_data = dims(2)*dims(3);
            n_convolution = n_wavelet+n_data-1;
            n_conv_pow2 = pow2(nextpow2(n_convolution));
            half_of_wavelet_size = (n_wavelet-1)/2;
            
            % For Laplacian
            X = [BV_Chanlocs_60.X]; Y = [BV_Chanlocs_60.Y]; Z = [BV_Chanlocs_60.Z];
            
            % Pick channel
            chans=[36,33,56];    % FCz, F5, F6
            
            for REFi=1:2
                if REFi==1, TAG='V';
                elseif REFi==2, TAG='L';
                    [EEG.data,~,~] = laplacian_perrinX(EEG.data,X,Y,Z,[],1e-6);
                end
                
                % Get FFT of data
                for chani=1:3
                    EEG_fft(chani,:) = fft(reshape(EEG.data(chans(chani),:,:),1,n_data),n_conv_pow2);
                end
                
                for fi=1:num_freqs
                    
                    wavelet = fft( exp(2*1i*pi*frex(fi).*t) .* exp(-t.^2./(2*(s(fi)^2))) , n_conv_pow2 );  % sqrt(1/(s(fi)*sqrt(pi))) *
                    
                    % convolution
                    for chani=1:3
                        temp_conv = ifft(wavelet.*EEG_fft(chani,:));
                        temp_conv = temp_conv(1:n_convolution);
                        temp_conv = temp_conv(half_of_wavelet_size+1:end-half_of_wavelet_size);
                        EEG_conv(chani,:,:) = reshape(temp_conv,dims(2),dims(3));
                        clear temp_conv;

                        % Common pre-EEG baseline
                        temp_BASE(chani,:) = mean(mean(abs(EEG_conv(chani,B1:B2,:)).^2,2),3);
                    end
                    
                    for idx=1:3
                        
                        if     idx==1,  idx_V=VECTOR(:,1)==201;  % STD
                        elseif idx==2,  idx_V=VECTOR(:,1)==200;  % TARG
                        elseif idx==3,  idx_V=VECTOR(:,1)==202;   % NOV
                        end
                        
                        for chani=1:3
                            temp_PWR = squeeze(mean(abs(EEG_conv(chani,T1:T2,idx_V)).^2,3));
                            
                            POWER(chani,fi,:,idx)   = 10* (  log10(temp_PWR') - log10(repmat(temp_BASE(chani,:),size(tx2disp,2),1))  );
                            ITPC(chani,fi,:,idx)    = abs(mean(exp(1i*(  angle(EEG_conv(chani,T1:T2,idx_V))  )),3));
                            if     chani==1, seed=1; targ=2;
                            elseif chani==2, seed=1; targ=3;
                            elseif chani==3, seed=2; targ=3;
                            end
                            ISPC(chani,fi,:,idx) = abs(mean(exp(1i*(  angle(EEG_conv(seed,T1:T2,idx_V)) - angle(EEG_conv(targ,T1:T2,idx_V))  )),3));
                            
                            clear temp_PWR;
                        end
                        clear idx_V ;
                    end
                    clear wavelet idx_V temp_BASE EEG_conv;
                end
                
                %%
                % $$$$$$$$$$$$$$$$$$$$$$$            $$$$$$$$$$$$$$$$$$$$$$$
                % $$$$$$$$$$$$$$$$$$$$$$$ Theta Topo
                % $$$$$$$$$$$$$$$$$$$$$$$            $$$$$$$$$$$$$$$$$$$$$$$
                
                topofrex=4.5;
                s=logspace(log10(3),log10(10),num_freqs)./(2*pi*topofrex);
                wavelet = fft( exp(2*1i*pi*frex(fi).*t) .* exp(-t.^2./(2*(s(fi)^2))) , n_conv_pow2 );  % sqrt(1/(s(fi)*sqrt(pi))) *
                
                seed=36;
                
                EEG_fft_4TOPO = fft(reshape(EEG.data(seed,:,:),1,n_data),n_conv_pow2);
                seed_EEG_conv_4TOPO = ifft(wavelet.*EEG_fft_4TOPO);
                seed_EEG_conv_4TOPO = seed_EEG_conv_4TOPO(1:n_convolution);
                seed_EEG_conv_4TOPO = seed_EEG_conv_4TOPO(half_of_wavelet_size+1:end-half_of_wavelet_size);
                seed_EEG_conv_4TOPO = reshape(seed_EEG_conv_4TOPO,dims(2),dims(3));
                clear EEG_fft_4TOPO ;
                
                % Common pre-EEG SEED baseline
                seed_BASE = mean(mean(abs(seed_EEG_conv_4TOPO(B1:B2,:)).^2,1),2);
                
                for chani=1:60
                    
                    EEG_fft_4TOPO = fft(reshape(EEG.data(chani,:,:),1,n_data),n_conv_pow2);
                    EEG_conv_4TOPO = ifft(wavelet.*EEG_fft_4TOPO);
                    EEG_conv_4TOPO = EEG_conv_4TOPO(1:n_convolution);
                    EEG_conv_4TOPO = EEG_conv_4TOPO(half_of_wavelet_size+1:end-half_of_wavelet_size);
                    EEG_conv_4TOPO = reshape(EEG_conv_4TOPO,dims(2),dims(3));

                    % Common pre-EEG baseline
                    temp_BASE = mean(mean(abs(EEG_conv_4TOPO(B1:B2,:)).^2,1),2);
                    
                    for idx=1:3
                        if     idx==1,  idx_V=VECTOR(:,1)==201;  % STD
                        elseif idx==2,  idx_V=VECTOR(:,1)==200;  % TARG
                        elseif idx==3,  idx_V=VECTOR(:,1)==202;   % NOV
                        end
                        
                        temp_PWR = squeeze(mean(abs(EEG_conv_4TOPO(T1:T2,idx_V)).^2,2));
                        POWER_TOPO(chani,:,idx) = 10* (  log10(temp_PWR) - log10(repmat(temp_BASE,size(tx2disp,2),1))  );
                        
                        S4cor=10* (  log10(abs(seed_EEG_conv_4TOPO(T1:T2,idx_V)).^2) -  log10(repmat(seed_BASE,size(tx2disp,2),sum(idx_V))) );
                        T4cor=10* (  log10(abs(EEG_conv_4TOPO(T1:T2,idx_V)).^2) -  log10(repmat(temp_BASE,size(tx2disp,2),sum(idx_V))) );
                        CORREL_TOPO(chani,:,idx)= diag(corr(S4cor',T4cor','type','Spearman'));
                        
                        SYNCH_TOPO(chani,:,idx) = abs(mean(exp(1i*(  angle(seed_EEG_conv_4TOPO(T1:T2,idx_V)) - angle(EEG_conv_4TOPO(T1:T2,idx_V))  )),2));
                        
                        clear idx_V temp_PWR  S4cor T4cor;
                    end
                    
                    clear EEG_fft_4TOPO  EEG_conv_4TOPO  TOPO_conv temp_BASE;
                end
                
                %%
                % $$$$$$$$$$$$$$$$$$$$$$$      $$$$$$$$$$$$$$$$$$$$$$$
                % $$$$$$$$$$$$$$$$$$$$$$$ ERPs
                % $$$$$$$$$$$$$$$$$$$$$$$      $$$$$$$$$$$$$$$$$$$$$$$
                
                % Filter
                dims=size(EEG.data);
                EEG.data=eegfilt(EEG.data,500,[],20);
                EEG.data=eegfiltfft(EEG.data,500,.1,[]);
                EEG.data=reshape(EEG.data,dims(1),dims(2),dims(3));
                
                % Basecor your ERPs here so they are pretty.
                EEG_BASE=squeeze( mean(EEG.data(:,find(tx==-200):find(tx==0),:),2) );
                for ai=1:dims(1)
                    EEG.data(ai,:,:)=squeeze(EEG.data(ai,:,:))-repmat( EEG_BASE(ai,:),dims(2),1 );
                end
                
                % Get ERPs
                for idx=1:3
                        if     idx==1,  idx_V=VECTOR(:,1)==201;  % STD
                        elseif idx==2,  idx_V=VECTOR(:,1)==200;  % TARG
                        elseif idx==3,  idx_V=VECTOR(:,1)==202;   % NOV
                        end
                    ERP(:,:,idx)=squeeze(mean(EEG.data(:,find(tx==-500):find(tx==1000),idx_V),3));
                    clear DATA_erp idx_V ;
                end
                
                save([savedir,num2str(subno),'_',num2str(sess),'_3AOB_TFandERPs_',TAG,'.mat'],'ERP','ISPC','ITPC','POWER','VECTOR','SYNCH_TOPO','TRL_ct','POWER_TOPO','CORREL_TOPO');
                
                clear ERP ISPC ITPC POWER RT;
            end
            
            clearvars -except datadir savedir FILEENDER BV_Chanlocs_60 D_DAT D_HDR D_ALL tx B1 B2 T1 T2 tx2disp si sess
        end
    end
end

%%
