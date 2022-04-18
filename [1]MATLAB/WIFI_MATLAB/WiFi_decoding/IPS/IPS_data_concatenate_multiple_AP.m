%%
%
%  To prepare the IPS training data
%  Creator: zzf
% CNN
% clear
% clc

%% Parameter Setting
fileDate = {'0417', '0421'};
targ = 'train';
% File Save Path
exportPath = '0423';

exportPath = strcat('~/ips/cnn_ips/data/',exportPath, '/DL_Data/');
tarData = strcat(targ, '_data');
tarLabel = strcat(targ, '_label');
dataPath = strcat(exportPath, targ, '/', tarData);
labelPath = strcat(exportPath, targ, '/', tarLabel);

numSample = 2000;
numTrain = 31;
numTest = 19;
numberAP = length(fileDate);

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
    fileDate_ = fileDate{nnn};
    importPath = strcat('~/ips/cnn_ips/data/', fileDate_, '/Matlab_Data/');
    importPath = strcat(importPath, targ, '/lltf/');
    
    %% Data import and concatenation
    comb = zeros(totalNum, 52);
    for i = 1:numPoint
        location = num2str(i);
        conc = strcat(importPath,location);
        load(conc)
        comb((i-1)*numSample+1:i*numSample,:) = csi_lltf;
    end
    
    %% Generate training data: Shuffle
    % Generate random number
    rng(23)
    ind_ran = randperm(totalNum);
    if isequal(targ, 'train')
        train_data(totalNum * (nnn-1)+1:totalNum * nnn,:) = comb(ind_ran',:);
    else
        test_data(totalNum * (nnn-1)+1:totalNum * nnn,:) = comb(ind_ran',:);
    end
    
    %% Generate Lable
    label = zeros(totalNum, 1);
    ss = 0;
    for k = 1 : numPoint
        label((k-1)*numSample+1:k*numSample,:) = ones(numSample,1) * ss;
        ss = ss + 1;
    end
    if isequal(targ, 'train')
        train_label(totalNum * (nnn-1)+1:totalNum * nnn) = label(ind_ran);
    else
        test_label(totalNum * (nnn-1)+1:totalNum * nnn) = label(ind_ran);
    end
end

%% Data Save
save(dataPath, tarData)
save(labelPath, tarLabel)

