%%
%
%  To prepare the IPS training data
%  Creator: zzf
% CNN
clear
% clc

%% Parameter Setting
% fileDate = {'0417', '0421'};
fileDate = {'0421'};

targ = 'test';
% File Save Path
exportPath = '0421';

exportPath = strcat('~/ips/cnn_ips/data/',exportPath, '/DL_Data/');
tarData = strcat(targ, '_data');
tarLabel = strcat(targ, '_label');
dataPath = strcat(exportPath, targ, '/', tarData);
labelPath = strcat(exportPath, targ, '/', tarLabel);

numSample = 850;
numTrain = 31;
numTest = 19;
numberAP = length(fileDate);

if isequal(targ, 'train')
    numPoint = numTrain;
    totalNum = numPoint*numSample;
    train_data = zeros(totalNum, 52);
    train_label = zeros(totalNum, 1);
else
    numPoint = numTest;
    totalNum = numPoint*numSample;
    test_data = zeros(totalNum, 52);
    test_label = zeros(totalNum, 1);
end

%% Data import and concatenation
for nnn = 1 : numberAP
    fileDate_ = fileDate{nnn};
    importPath = strcat('~/ips/cnn_ips/data/', fileDate_, '/Matlab_Data/');
    importPath = strcat(importPath, targ, '/lltf/');
    
    comb = zeros(totalNum, 52);
    for i = 1:numPoint
        location = num2str(i);
        conc = strcat(importPath,location);
        load(conc)
        comb((i-1)*numSample+1:i*numSample,:) = csi_lltf(1:numSample, :);
    end
    % Adding up two data sets to eliminate symmetry
    if isequal(targ, 'train')
        train_data= train_data + comb;
    else
        test_data = test_data + comb;
    end
end

%% Generate training data: Shuffle
% Generate random number
rng(23)
ind_ran = randperm(totalNum);
if isequal(targ, 'train')
    train_data= train_data(ind_ran',:);
else
    test_data = test_data(ind_ran',:);
end

%% Generate Lable
label = zeros(totalNum, 1);
ss = 0;
for k = 1 : numPoint
    label((k-1)*numSample+1:k*numSample,:) = ones(numSample,1) * ss;
    ss = ss + 1;
end

if isequal(targ, 'train')
    train_label= label(ind_ran);
else
    test_label = label(ind_ran);
end

%% Data Save
save(dataPath, tarData)
save(labelPath, tarLabel)

