%%
%
%  PREFORMATTED
%  TEXT
%
% clear
% clc

%% Load predicted label
dataType = 'test';
rootPath = 'C:\Users\HYPC300\Desktop\ips_data\';
dataInput = '0421_31pts_ap2';
dataOutput = '0506_timeDomain_64sample_ap2';
sampleFileName = "\20M_wifi.txt";

numTrainPoint = 31;
numTestPoint = 19;
numFrame = 2000;

%% Path setting
inputPath = strcat(rootPath, dataInput, '\', dataType, '\');
outputDataPath = strcat(rootPath, dataOutput, '\', dataType, '\data');
outputFigurePath = strcat(rootPath, dataOutput, '\', dataType, '\figure');

if isequal (dataType, 'train')
    numData = numTrainPoint; % numTest
else
    numData = numTestPoint; % numTest
end

%% WiFi parameter initialization
cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS', 0);

if cfgVHT.MCS ==8
    sampleNum = 2960;
else
    sampleNum = 960;
end
%% VHT configuration
ind_Data = wlanFieldIndices(cfgVHT,'VHT-Data');
ind_VHTltf = wlanFieldIndices(cfgVHT,'VHT-LTF');
ind_Lltf = wlanFieldIndices(cfgVHT,'L-LTF');

%%
numOutlier = zeros(numData, 1);
epss = zeros(numData, 1);

for indexNum = 1 : numData
    close all;
    % Load the data from txt file for each point
    Input = IPS_dataLoad(inputPath, indexNum, sampleFileName);
    %%  Synchronization
    % The number of synched frames
    %     Synchronization, and return the corresponding frame in the packet
    rxSig_collect = double(IPS_synchronization(Input, cfgVHT, numFrame));
    
    %%  Plot the time domain waveform of L-LTF (160samples)
    lltf_time = abs(rxSig_collect(ind_Lltf(1):ind_Lltf(2),:));
    lltf_time_reCP = lltf_time(33:end, :, :);
    
    lltf_time_rs = reshape(lltf_time_reCP, 64, 2, numFrame);
    lltf_time_ave = (lltf_time_rs(:,1,:) + lltf_time_rs(:,2,:))/2;
    lltf_time_ave = squeeze(lltf_time_ave);
    
    % Conform to DL data
    csi_lltf = lltf_time_ave';
    
    % Plot the channel using VHT-LTF
    figure(indexNum)
    hold on
    for k = 1: numFrame
        plot(lltf_time_ave(:,k));
    end
    
    % Save lltf variable to .mat file
    strr = num2str(indexNum);
    path_ = strcat(outputDataPath, '\lltf\', strr);
    save(path_, 'csi_lltf') ;
    
    % Save lltf figure
    path_ = strcat(outputFigurePath, '\lltf\', strr);
    saveas(indexNum, path_, 'png' )
end