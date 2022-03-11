%%
%
%  PREFORMATTED
%  TEXT
%
clear

%% Load predicted label
rootPath = 'C:\Users\HYPC300\Desktop\ips_data\';
dataInput = '0508_device_based_ap4';
dataOutput = '0508_device_based_ap4';
numTrainPoint = 5;
numTestPoint = 4;
numFrame = 2000;
numDiscard = 500;

%% Path setting
dataType = 'test';
inputPath = strcat(rootPath, dataInput, '\', dataType, '\data\lltf\');
outputDataPath = strcat(rootPath, dataOutput, '\', dataType, '\data\lltf\');
outputFigurePath = strcat(rootPath, dataOutput, '\', dataType, '\figure\lltf\');

if isequal (dataType, 'train')
    numData = numTrainPoint; % numTest
else
    numData = numTestPoint; % numTest
end

%%
numOutlier = zeros(numData, 1);
epss = zeros(numData, 1);

for i = 1 : numData
    close all;
    % the parameter for dbscan
    eps = 50;
    numm = num2str(i);
    fileName = strcat(inputPath, numm);
    load(fileName);
    % apply dbscan algorithm
    chEst_lltf_ind = dbscan(csi_lltf, eps, 30);
    discardd = sum(chEst_lltf_ind ~= 1);
    while true
        if discardd >= numDiscard
            eps = eps + 1;
            chEst_lltf_ind = dbscan(csi_lltf, eps, 30);
            discardd = sum(chEst_lltf_ind ~= 1);
        else
            break;
        end
    end
    % recorder the number of discarded samples
    epss(i) = eps;
    numOutlier(i) = discardd;
    
    figure(i)
    hold on
    indd = find(chEst_lltf_ind == 1);
    chEst_lltf_abs_den = csi_lltf(indd, :);
    for k = 1: length(indd)
        plot(chEst_lltf_abs_den(k, :));
    end
    csi_lltf = chEst_lltf_abs_den;
    
    % Save lltf variable to .mat file
    path_ = strcat(outputDataPath, numm);
    save(path_, 'csi_lltf') ;
    
    % Save lltf figure
    path_ = strcat(outputFigurePath, numm);
    saveas(i, path_, 'png' )
end