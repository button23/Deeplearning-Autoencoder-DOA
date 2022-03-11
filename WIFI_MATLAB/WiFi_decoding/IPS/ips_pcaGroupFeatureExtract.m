%%
%
%  PREFORMATTED
%  TEXT
%
clear;
close all
clc
%% Load predicted label
rootPath = 'C:\Users\HYPC300\Desktop\ips_data\';
dataInput = '0421_wifi_ips_data_denoise_ap2';
dataOutput = '0430_wifi_ips_denoisedata_pca_ap2';
numTestPoint = 19;
numFrame = 1500;
pcaRatio = 0.99;
numSCs = 52;

targ = 'test';

%% Path setting
inputPath = strcat(rootPath, dataInput, '\', targ, '\data\lltf\');
outputDataPath = strcat(rootPath, dataOutput, '\', targ, '\');
pcagPath = strcat(rootPath, 'pcag');

%% Data loading
if isequal(targ, 'train')
    k_record = zeros(numTestPoint, 1);
    % load pca group (binding training points to extract mutual features)
    load(pcagPath)
    % load training data
    for i = 1 : numTestPoint
        % inside pcag defines the specific index of training points for every
        % training point
        subgroup = pcag{i};
        concCSI = zeros(numFrame, numSCs * length(subgroup));
        % load neighboring training data for each testing point
        for j = 1 : length(subgroup)
            numm = num2str(subgroup(j));
            fileName = strcat(inputPath, numm);
            load(fileName);
            % Concatenate the neighboring training samples
            concCSI(:, numSCs*(j-1)+1:numSCs*j) = csi_lltf(1:numFrame,:);
        end
        % apply pca algorithm on concatenated samples
        %% Apply pca algorithm
        [signals, PC, V] = pca2(concCSI.');
        csi_lltf = signals(1:numSCs,:)';
        
        % Save lltf variable to .mat file
        nnn = num2str(i);
        path_ = strcat(outputDataPath, nnn);
        save(path_, 'csi_lltf') ;
    end
else
    % load testing data
    for i = 1 : numTestPoint
        numm = num2str(i);
        fileName = strcat(inputPath, numm);
        load(fileName);
        %% Apply pca algorithm
        [signals, PC, V] = pca2(csi_lltf.');
        csi_lltf = signals(1:numSCs,:)';
        
        % Save lltf variable to .mat file
        nnn = num2str(i);
        path_ = strcat(outputDataPath, nnn);
        save(path_, 'csi_lltf') ;
    end
end
%%