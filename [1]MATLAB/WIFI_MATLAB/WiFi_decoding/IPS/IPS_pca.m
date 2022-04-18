%%
%
%  To prepare the IPS training data
%  Creator: zzf
% CNN
clear
clc

%% Parameter Setting
fileDate = {'0421_wifi_ips_data_denoise_ap2'}; % 0421_wifi_ips_data_denoise_ap2   0421_wifi_ips_data_ap2
targ = 'train';
% File Save Path
exportPath = '0421_wifi_ips_data_ap2';

exportPath = strcat('C:/Users/HYPC300/Desktop/ips_data/',exportPath);
tarData = strcat(targ, '_data');
tarLabel = strcat(targ, '_label');
dataPath = strcat(exportPath, targ, '/data_pca/', tarData);
labelPath = strcat(exportPath, targ, '/data_pca/', tarLabel);

numSample = 1500;
numTrain = 31;
numTest = 19;
numberAP = length(fileDate);
pcaRatio = 0.999;  % decides the degree of reduction on dimentions

if isequal(targ, 'train')
    numPoint = numTrain;
    totalNum = numPoint*numSample;
    train_data = zeros(totalNum * numberAP, 52);
    train_label = zeros(totalNum * numberAP, 1);
else
    numPoint = numTest;
    totalNum = numPoint*numSample;
    test_data = zeros(totalNum * numberAP, 52);
    test_label = zeros(totalNum * numberAP, 1);
end

%% Load Path
for nnn = 1 : numberAP
    k_record = zeros(numPoint, 1);
    fileDate_ = fileDate{nnn};
    importPath = strcat('C:/Users/HYPC300/Desktop/ips_data/', fileDate_, '/', targ, '/data/lltf/');
    
    %% Data import and concatenation
    comb = zeros(totalNum, 52);
    for i = 1:numPoint
        location = num2str(i);
        conc = strcat(importPath, location);
        load(conc)
        csi_lltf = csi_lltf(1:numSample, :);
        %% Apply pca algorithm
%         [signals, PC, V] = pca2(csi_lltf.');
        [signals, PC, V] = pca2(c1.');
        k = 1;
        flag = 1;
        while flag==1
            rati = (sum(V(1:k))/sum(V));
            if rati >= pcaRatio
                flag = 0;
            else
                k = k+1;
            end
        end
        k_record(i) = k;
%         comb((i-1)*numSample+1:i*numSample,:) = csi_lltf;
    end
    
%     %% Generate training data: Shuffle
%     % Generate random number
%     rng(23)
%     ind_ran = randperm(totalNum);
%     if isequal(targ, 'train')
%         train_data(totalNum * (nnn-1)+1:totalNum * nnn,:) = comb(ind_ran',:);
%     else
%         test_data(totalNum * (nnn-1)+1:totalNum * nnn,:) = comb(ind_ran',:);
%     end
%     
%     %% Generate Lable
%     label = zeros(totalNum, 1);
%     ss = 0;
%     for k = 1 : numPoint
%         label((k-1)*numSample+1:k*numSample,:) = ones(numSample,1) * ss;
%         ss = ss + 1;
%     end
%     if isequal(targ, 'train')
%         train_label(totalNum * (nnn-1)+1:totalNum * nnn) = label(ind_ran);
%     else
%         test_label(totalNum * (nnn-1)+1:totalNum * nnn) = label(ind_ran);
%     end
end

%% Data Save
% save(dataPath, tarData)
% save(labelPath, tarLabel)
