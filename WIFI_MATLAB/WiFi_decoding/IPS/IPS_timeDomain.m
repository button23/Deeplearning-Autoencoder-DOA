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
dataOutput = '0501_timedomain_ap2';
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

%%
numOutlier = zeros(numData, 1);
epss = zeros(numData, 1);

for indexNum = 1 : numData
    close all;
    % Load the data from txt file for each point
    Input = IPS_dataLoad(inputPath, indexNum, sampleFileName);
    %%  Synchronization
    % The number of synched frames
    rxSig_collect = zeros(sampleNum, numFrame);
    %     numDec = 2; % Decoding the nth frame from total frame
    %     Synchronization, and return the corresponding frame in the packet
    rxSig_collect(: , 1 : numFrame) = IPS_synchronization(Input, cfgVHT, numFrame);
    
    %% VHT configuration
    ind_Data = wlanFieldIndices(cfgVHT,'VHT-Data');
    ind_VHTltf = wlanFieldIndices(cfgVHT,'VHT-LTF');
    ind_Lltf = wlanFieldIndices(cfgVHT,'L-LTF');
    
    %%
    lltf_time = abs(rxSig_collect(ind_Lltf(1):ind_Lltf(2),:));
    vhtltf_time = abs(rxSig_collect(ind_VHTltf(1):ind_VHTltf(2),:));
    
    %%
    % Plot the channel using L-LTF
    figure(indexNum)
    hold on
    for k = 1: numFrame
        plot(lltf_time(:,k));
    end
    % Conform to DL data
    csi_lltf = lltf_time';
    
    % Save lltf variable to .mat file
    strr = num2str(indexNum);
    path_ = strcat(outputDataPath, '\lltf\', strr);
    save(path_, 'csi_lltf') ;
    
    % Save lltf figure
    path_ = strcat(outputFigurePath, '\lltf\', strr);
    saveas(indexNum, path_, 'png' )
    
    %%
    % Plot the channel using VHT-LTF
    figure(indexNum+100)
    hold on
    for k = 1: numFrame
        plot(vhtltf_time(:,k));
    end
    % Conform to DL data
    csl_vhtltf = vhtltf_time';
    
    % Save vhtltf variable to .mat file
    path_ = strcat(outputDataPath, '\vhtltf\', strr);
    save(path_, 'csl_vhtltf') ;
    
    % Save vhtltf figure
    path_ = strcat(outputFigurePath, '\vhtltf\', strr);
    saveas(indexNum+100, path_, 'png' )
end