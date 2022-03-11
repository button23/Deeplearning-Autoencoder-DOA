%% Initialization
clear;
clc;
%% WiFi parameter initialization
cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS', 0);

if cfgVHT.MCS ==8
    sampleNum = 2960;
else
    sampleNum = 960;
end

%% Load RF Raw data
trainPath = 'C:\Users\HYPC300\Desktop\ips_data\0508_device_based_ap4\train\';
testPath = 'C:\Users\HYPC300\Desktop\ips_data\0508_device_based_ap4\test\';
% trainPath = 'C:\Users\HYPC300\Desktop\ips_data\not_human\';
% sampleFileName = ["\a1.txt","\a2.txt","\a3.txt","\a4.txt","\a5.txt"];
% sampleFileName = ["\20M_wifi.txt"];
sampleFileName = ["\a1.txt"];

pathFile = testPath;

numTrainPoint = 5;
numTestPoint = 4 ;
numFile = 1; %5
numFrame = 2000;

if isequal(pathFile, trainPath)
    numPoint = numTrainPoint;
else
    numPoint = numTestPoint;
end

% Specify the index of the point (FROM 62 training point, or 43 test point)
% indexNum = 1;

for indexNum = 1 : numPoint
    close all;
    % Load the data from txt file for each point
    Input = IPS_dataLoad(pathFile, indexNum, sampleFileName);
    
    %%  Synchronization
    % The number of synched frames
    rxSig_collect = zeros(sampleNum, numFrame);
    %     numDec = 2; % Decoding the nth frame from total frame
    %     Synchronization, and return the corresponding frame in the packet
        rxSig_collect(: ,1 : numFrame) = IPS_synchronization(Input, cfgVHT, numFrame);
    
    %% VHT configuration
    ind_Data = wlanFieldIndices(cfgVHT,'VHT-Data');
    ind_VHTltf = wlanFieldIndices(cfgVHT,'VHT-LTF');
    ind_Lltf = wlanFieldIndices(cfgVHT,'L-LTF');
    
    %% Data recover (WLAN function)  ind_VHTltf
    demodulateLLTF = wlanLLTFDemodulate(rxSig_collect(ind_Lltf(1):ind_Lltf(2),:), cfgVHT);
    chEst_lltf = wlanLLTFChannelEstimate(demodulateLLTF,cfgVHT);
    ltfDemodSig = wlanVHTLTFDemodulate(rxSig_collect(ind_VHTltf(1):ind_VHTltf(2),:), cfgVHT);
%     chEst_vhtltf = wlanVHTLTFChannelEstimate(ltfDemodSig,cfgVHT);
    chEst_lltf_s = squeeze(chEst_lltf);
%     chEst_vhtltf_s = squeeze(chEst_vhtltf);
    
    % Absolute value
    chEst_lltf_abs = abs(chEst_lltf_s);
%     chEst_vhtltf_abs = abs(chEst_vhtltf_s);
    
    % Conform to DL data
    csi_lltf = abs(chEst_lltf_abs.');
%     csl_vhtltf = abs(chEst_vhtltf_abs.');
    
    % Plot the channel using L-LTF
    figure(indexNum)
    hold on
    for k = 1: numFrame
        plot(chEst_lltf_abs(:,k));
    end
    
    % Save lltf variable to .mat file
    strr = num2str(indexNum);
    path_ = strcat(pathFile, 'data\lltf\', strr);
    save(path_, 'csi_lltf') ;
    
    % Save lltf figure
    path_ = strcat(pathFile, 'figure\lltf\', strr);
    saveas(indexNum, path_, 'png' )
    
%     % Plot the channel using VHT-LTF
%     figure(indexNum+100)
%     hold on
%     for k = 1: numFrame
%         plot(chEst_vhtltf_abs(:,k));
%     end
%     
%     % Save vhtltf variable to .mat file
%     path_ = strcat(pathFile, 'data\vhtltf\', strr);
%     save(path_, 'csl_vhtltf') ;
%     
%     % Save vhtltf figure
%     path_ = strcat(pathFile, 'figure\vhtltf\', strr);
%     saveas(indexNum+100, path_, 'png' )
end
