%%
% 
%  To prepare the IPS training data
%  Creator: zzf
% CNN
% clear
% clc
%% Load Path and Save Path
targ = 'train';

tarData = strcat(targ, '_data');
tarLabel = strcat(targ, '_label');

importPath = '~/ips/cnn_ips/data/0417/Matlab_Data/';
importPath = strcat(importPath, targ, '/data/lltf/');
dataPath = '~/ips/cnn_ips/data/0417/DL_Data/';
dataPath = strcat(dataPath, targ, '/', tarData);
labelPath = '~/ips/cnn_ips/data/0417/DL_Data/';
labelPath = strcat(labelPath, targ, '/', tarLabel);

%% Parameter Setting
numSample = 2000;
numTrain = 31;
numTest = 19;

if isequal(targ, 'train')
    numPoint = numTrain;
else
    numPoint = numTest;
end

%% Data import and concatenation
totalNum = numPoint*numSample;
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
    train_data = comb(ind_ran',:);
else
    test_data = comb(ind_ran',:);
end
    save(dataPath, tarData)

%% Generate Lable
label = zeros(totalNum, 1);
ss = 0;

for k = 1 : numPoint
    label((k-1)*numSample+1:k*numSample,:) = ones(numSample,1) * ss;
    ss = ss + 1;
end
if isequal(targ, 'train')
    train_label = label(ind_ran);
else
    test_label = label(ind_ran);
end
save(labelPath, tarLabel)
