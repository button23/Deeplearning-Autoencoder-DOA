close all;
clear;
clc;

numTrain = 2000;
numTest = numTrain * 0.2;
validLength = 560;

snrr_record_train = zeros(numTrain,3);
snrr_record_test = zeros(numTest,3);

train_data_ = zeros(numTrain*3, validLength);
test_data_ = zeros(numTest*3, validLength);

%% wlan waveform
cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS',3); % 16qam 1/2 data rate
cfgHT = wlanHTConfig('ChannelBandwidth','CBW20','MCS',3); % 16qam 1/2 data rate
cfgnonHT = wlanNonHTConfig('ChannelBandwidth','CBW20','MCS',4); % 16qam 1/2 data rate
cfgVHT.APEPLength = 150;
cfgHT.PSDULength = 150;
cfgnonHT.PSDULength = 200;


% Random snr
% snrr = 5;
% snrr = randi([-5,20], (numTrain+numTest),1);
% for i=1:length(mcss)
    tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-F');


% pfOffset = comm.PhaseFrequencyOffset('SampleRate',20e6,'FrequencyOffsetSource','Input port');

for mjk =1:(numTrain + numTest)
    % Data generation
    PSDU_VHT = randi([0 1],cfgVHT.PSDULength*8,1);
    PSDU_HT = randi([0 1],cfgVHT.PSDULength*8,1);
    PSDU_nonHT = randi([0 1],cfgVHT.PSDULength*8,1);
    
    wholeWaveform_VHT = wlanWaveformGenerator(PSDU_VHT,cfgVHT);
    wholeWaveform_HT = wlanWaveformGenerator(PSDU_HT,cfgHT);
    wholeWaveform_nonHT = wlanWaveformGenerator(PSDU_nonHT,cfgnonHT);
    
%     rng(2020)

    %% Data Convertion
    % Create random bit stream
    %       rng(25)
    % Channel and AWGN
    snr_vht = randi([0 30]);
    chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR', snr_vht);
    rxSig_VHT = chNoise(tgacChan(wholeWaveform_VHT));  %% wholeWaveform_scaling  wholeWaveform
%         rxSig_VHT = chNoise(wholeWaveform_VHT);  %% wholeWaveform_scaling  wholeWaveform

    
    snr_ht = randi([0 30]);
    chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR', snr_ht);
    rxSig_HT = chNoise(tgacChan(wholeWaveform_HT));
%         rxSig_HT = chNoise(wholeWaveform_HT);

    snr_nonht = randi([0 30]);
    chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR', snr_nonht);
    rxSig_nonHT = chNoise(tgacChan(wholeWaveform_nonHT));
%         rxSig_nonHT = chNoise(wholeWaveform_nonHT);
    
    if mjk <= numTrain
        train_data_(mjk,:) = rxSig_VHT(1:validLength);
        snrr_record_train(mjk,1) = snr_vht;
        
        train_data_((numTrain+mjk),:) = rxSig_HT(1:validLength);
        snrr_record_train(mjk,2) = snr_ht;
        
        train_data_((numTrain * 2 + mjk),:) = rxSig_nonHT(1:validLength);
        snrr_record_train(mjk,3) = snr_nonht;
    else
        test_data_(mjk-numTrain,:) = rxSig_VHT(1:validLength);
        % Record the SNR of the test data.
        snrr_record_test(mjk-numTrain, 1) = snr_vht;
        test_data_((mjk-numTrain+numTest),:) = rxSig_HT(1:validLength);
        % Record the SNR of the test data.
        snrr_record_test(mjk-numTrain, 2) = snr_ht;
        test_data_((mjk-numTrain+numTest*2),:) = rxSig_nonHT(1:validLength);
        % Record the SNR of the test data.
        snrr_record_test(mjk-numTrain, 3) = snr_nonht;
    end
end

% Training set label
train_label_ = [zeros(numTrain,1); ones(numTrain,1); 2*ones(numTrain,1)];
% Test set label
test_label_ = [zeros(numTest,1); ones(numTest,1); 2*ones(numTest,1)];


%% Shuffle training data
snrr_record_train_ = reshape(snrr_record_train, numel(snrr_record_train),1);
shuf_ = [train_data_,train_label_,snrr_record_train_];
indd_ = randperm(3*numTrain);
test_datass_ = shuf_(indd_,:);
train_data = test_datass_(:,1:validLength);
train_label = test_datass_(:,end-1);
train_snr = test_datass_(:,end);

% test_data = train_data;
% test_label = train_label;
% test_snr = train_snr;


%% Shuffle test data
snrr_record = reshape(snrr_record_test,numel(snrr_record_test),1);
shuf = [test_data_,test_label_,snrr_record];
indd = randperm(3*numTest);
test_datass = shuf(indd,:);
test_data = test_datass(:,1:validLength);
test_label = test_datass(:,end-1);
test_snr = test_datass(:,end);





