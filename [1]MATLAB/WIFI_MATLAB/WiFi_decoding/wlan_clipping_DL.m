close all;
clear;
clc;

numTrain = 1000;
numTest = numTrain * 0;
validLength = 640;

SNRR = -10:20;
snrr_record_train1 = zeros(numTrain * length(SNRR),1);
snrr_record_train2 = zeros(numTrain * length(SNRR),1);
snrr_record_train3 = zeros(numTrain * length(SNRR),1);

snr_drop = zeros(numTrain,3);

snrr_record_test1 = zeros(numTest * length(SNRR),1);
snrr_record_test2 = zeros(numTest * length(SNRR),1);
snrr_record_test3 = zeros(numTest * length(SNRR),1);

train_data_1 = zeros(numTrain * length(SNRR), validLength);
train_data_2 = zeros(numTrain * length(SNRR), validLength);
train_data_3 = zeros(numTrain * length(SNRR), validLength);

test_data_1= zeros(numTest * length(SNRR), validLength);
test_data_2= zeros(numTest * length(SNRR), validLength);
test_data_3= zeros(numTest * length(SNRR), validLength);

%% wlan waveform
cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS',0); % bpsk 1/2 data rate 3
cfgHT = wlanHTConfig('ChannelBandwidth','CBW20','MCS',0); % 16qam 1/2 data rate 3
cfgnonHT = wlanNonHTConfig('ChannelBandwidth','CBW20','MCS',0); % 16qam 1/2 data rate 4
cfgVHT.APEPLength = 150;
cfgHT.PSDULength = 150;
cfgnonHT.PSDULength = 200;


% Random snr
% snrr = 5;
% snrr = randi([-5,20], (numTrain+numTest),1);
% for i=1:length(mcss)
kk1 = 0;
kk2 = 0;
kk3 = 0;
threshold = 0.5;
for snrrr = 1:length(SNRR)
    
    
    if snrrr == 10
        fprintf('Progress 30\n');
    elseif snrrr == 19
        fprintf('Progress 60\n');
    elseif snrrr == 28
        fprintf('Progress 90\n');
    end
    
    % pfOffset = comm.PhaseFrequencyOffset('SampleRate',20e6,'FrequencyOffsetSource','Input port');
    %% Iteration
    for mjk =1:(numTrain + numTest)
        
%                 tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-E','LargeScaleFadingEffect','Pathloss and shadowing');
        tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-B');
        
        % Data generation
        PSDU_VHT = randi([0 1],cfgVHT.PSDULength*8,1);
        wholeWaveform_VHT = wlanWaveformGenerator(PSDU_VHT,cfgVHT);
        
        
        % clipping
         clipping_ratio = 1.0;
         clipped_VHT = wholeWaveform_VHT;
         clip_idx_VHT = abs(clipped_VHT) > clipping_ratio;
         
         for i=1:length(clipped_VHT)
             if clip_idx_VHT(i)
                 clipped_VHT(i) = (clipped_VHT(i)*clipping_ratio)./abs(wholeWaveform_VHT(i));
             end
         end
         
        
        %% Data Convertion
        % Create random bit stream
        %       rng(25)
        % Channel and AWGN
        %% VHT part
        %         snr_vht = randi([-10 20]);
        chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
            'SNR', SNRR(snrrr));
        rxSig_VHT = chNoise(tgacChan(clipped_VHT));  %% wholeWaveform_scaling  wholeWaveform
        %                 rxSig_VHT = chNoise((wholeWaveform_VHT));  %% wholeWaveform_scaling  wholeWaveform
                
%         %Syncronization
%         [startOffset_vht,M] = wlanSymbolTimingEstimate(rxSig_VHT, ...
%             cfgVHT.ChannelBandwidth,threshold);
%             startOffset_vht_rec(mjk) = startOffset_vht;
%         if startOffset_vht < 0 || startOffset_vht > 100
%             kk1 = kk1+1;
%             if mjk <= numTrain
%                 train_data_1((mjk+numTrain*(snrrr-1)),:) = 100;
%                 snrr_record_train1((mjk+numTrain*(snrrr-1))) = 100;
%             else
%                 test_data_1((mjk+numTrain*(snrrr-1))-numTrain,:) = 100;
%                 snrr_record_test1((mjk+numTrain*(snrrr-1))-numTrain) = 100;
%             end
%             continue
%         end
%         
%         rxSig_VHT_sync = rxSig_VHT(startOffset_vht+1:startOffset_vht+560);
        rxSig_VHT_sync = rxSig_VHT;
        if mjk <= numTrain
            train_data_1((mjk+numTrain*(snrrr-1)),:) = rxSig_VHT_sync(1:validLength);
            snrr_record_train1((mjk+numTrain*(snrrr-1))) = SNRR(snrrr);
        else
            test_data_1((mjk+numTrain*(snrrr-1))-numTrain,:) = rxSig_VHT_sync(1:validLength);
            snrr_record_test1((mjk+numTrain*(snrrr-1))-numTrain) = SNRR(snrrr);
        end
        
    end
    
    for mjk =1:(numTrain + numTest)
%                 tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-A','LargeScaleFadingEffect','Pathloss and shadowing');
        tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-B');
        
        PSDU_HT = randi([0 1],cfgVHT.PSDULength*8,1);
        wholeWaveform_HT = wlanWaveformGenerator(PSDU_HT,cfgHT);
        
        % clipping
         clipping_ratio = 1.0;
         clipped_HT = wholeWaveform_HT;
         clip_idx_HT = abs(clipped_HT) > clipping_ratio;
         
         for i=1:length(clipped_HT)
             if clip_idx_HT(i)
                 clipped_HT(i) = (clipped_HT(i)*clipping_ratio)./abs(wholeWaveform_HT(i));
             end
         end
        
        %% HT part
        %         snr_ht = randi([-10 20]);
        
        chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
            'SNR', SNRR(snrrr));
        rxSig_HT = chNoise(tgacChan(clipped_HT));
        
%         %Syncronization
%         startOffset_ht = wlanSymbolTimingEstimate(rxSig_HT, ...
%             cfgHT.ChannelBandwidth,threshold);
%         if startOffset_ht < 0 || startOffset_ht > 100
%             kk2 = kk2+1;
%             if mjk <= numTrain
%                 train_data_2((mjk+numTrain*(snrrr-1)),:) = 100;
%                 snrr_record_train2((mjk+numTrain*(snrrr-1))) = 100;
%             else
%                 test_data_2((mjk+numTrain*(snrrr-1))-numTrain,:) = 100;
%                 snrr_record_test2((mjk+numTrain*(snrrr-1))-numTrain) = 100;
%             end
%             continue
%         end
%         
%         rxSig_HT_sync = rxSig_HT(startOffset_ht+1:startOffset_ht+560);
        rxSig_HT_sync = rxSig_HT;
        if mjk <= numTrain
            train_data_2((mjk+numTrain*(snrrr-1)),:) = rxSig_HT_sync(1:validLength);
            snrr_record_train2((mjk+numTrain*(snrrr-1))) = SNRR(snrrr);
        else
            test_data_2((mjk+numTrain*(snrrr-1))-numTrain,:) = rxSig_HT_sync(1:validLength);
            snrr_record_test2((mjk+numTrain*(snrrr-1))-numTrain) = SNRR(snrrr);
        end
    end
    
    for mjk =1:(numTrain + numTest)
%                 tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-A','LargeScaleFadingEffect','Pathloss and shadowing');
        tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-B');
        
        PSDU_nonHT = randi([0 1],cfgVHT.PSDULength*8,1);
        wholeWaveform_nonHT = wlanWaveformGenerator(PSDU_nonHT,cfgnonHT);
        
        % clipping
         clipping_ratio = 1.0;
         clipped_nonHT = wholeWaveform_nonHT;
         clip_idx_nonHT = abs(clipped_nonHT) > clipping_ratio;
         
         for i=1:length(clipped_nonHT)
             if clip_idx_nonHT(i)
                 clipped_nonHT(i) = (clipped_nonHT(i)*clipping_ratio)./abs(wholeWaveform_nonHT(i));
             end
         end
        
        %% Non-HT part
        %         snr_nonht = randi([-10 20]);
        chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
            'SNR', SNRR(snrrr));
        rxSig_nonHT = chNoise(tgacChan(clipped_nonHT));
        
%         % Syncronization
%         startOffset_nonht = wlanSymbolTimingEstimate(rxSig_nonHT, ...
%             cfgnonHT.ChannelBandwidth,threshold);
%         
%         if startOffset_nonht < 0 || startOffset_nonht > 100
%             kk3 = kk3+1;
%             
%             if mjk <= numTrain
%                 train_data_3((mjk+numTrain*(snrrr-1)),:) = 100;
%                 snrr_record_train3((mjk+numTrain*(snrrr-1))) = 100;
%             else
%                 test_data_3((mjk+numTrain*(snrrr-1))-numTrain,:) = 100;
%                 snrr_record_test3((mjk+numTrain*(snrrr-1))-numTrain) = 100;
%             end
%             continue
%         end
%         
%         rxSig_nonHT_sync = rxSig_nonHT(startOffset_nonht+1:startOffset_nonht+560);
        rxSig_nonHT_sync = rxSig_nonHT;
        if mjk <= numTrain
            train_data_3((mjk+numTrain*(snrrr-1)),:) = rxSig_nonHT_sync(1:validLength);
            snrr_record_train3((mjk+numTrain*(snrrr-1))) = SNRR(snrrr);
        else
            test_data_3((mjk+numTrain*(snrrr-1))-numTrain,:) = rxSig_nonHT_sync(1:validLength);
            snrr_record_test3((mjk+numTrain*(snrrr-1))-numTrain) = SNRR(snrrr);
        end
        
    end
end


%% Training set label
train_label_ = [zeros(length(train_data_1(:,1)),1); ones(length(train_data_2(:,1)),1); 2*ones(length(train_data_3(:,1)),1)];
% Test set label
test_label_ = [zeros(length(test_data_1(:,1)),1); ones(length(test_data_1(:,1)),1); 2*ones(length(test_data_1(:,1)),1)];


%% Shuffle training data
% snrr_record_train_ = reshape(snrr_record_train, numel(snrr_record_train),1);
train_data_ = [train_data_1;train_data_2;train_data_3];
snrr_record_train_ = [snrr_record_train1;snrr_record_train2;snrr_record_train3];

firstCol = train_data_(:,1);
index_train = find(firstCol == 100);

train_data_(index_train,:)=[];
snrr_record_train_(index_train,:) = [];
train_label_(index_train,:) = [];

shuf_ = [train_data_,train_label_,snrr_record_train_];
indd_ = randperm(length(shuf_(:,1)));
test_datass_ = shuf_(indd_,:);

train_data = test_datass_(:,1:validLength);
train_label = test_datass_(:,end-1);
train_snr = test_datass_(:,end);

test_data = train_data;
test_label = train_label;
test_snr = train_snr;


%% Shuffle test data
% test_data_ = [test_data_1;test_data_2;test_data_3];
% snrr_record = [snrr_record_test1;snrr_record_test2;snrr_record_test3];
% 
% firstCol_ = test_data_(:,1);
% index_test = find(firstCol_ == 100);
% 
% test_data_(index_test,:)=[];
% snrr_record(index_test,:) = [];
% test_label_(index_test,:) = [];
% 
% shuf = [test_data_,test_label_,snrr_record];
% indd = randperm(length(shuf(:,1)));
% test_datass = shuf(indd,:);
% test_data = test_datass(:,1:validLength);
% test_label = test_datass(:,end-1);
% test_snr = test_datass(:,end);
