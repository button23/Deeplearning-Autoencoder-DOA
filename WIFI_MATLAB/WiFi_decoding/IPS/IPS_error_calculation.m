%%
%
%  PREFORMATTEDconfusedPosition
%  TEXT
%
clear
close all
%% Load predicted labelconfusedPosition
rootPath = '~/ips/cnn_ips/data/0421';
importPreLabel = strcat(rootPath, '/predict_data/predict_label.mat');
importTestLabel = strcat(rootPath,'/DL_Data/test/test_label.mat');
importErrorTable = strcat(rootPath,'/Matlab_Data/err_table.mat');
load(importPreLabel)
load(importTestLabel)
load(importErrorTable)

nClass = 19;
nSample = 1000;
%% Save path
resultMeanSavePath = strcat(rootPath,'/result/mean_error.mat');
resultStdSavePath = strcat(rootPath,'/result/std_error.mat');
resultMeanPPointSavePath = strcat(rootPath,'/result/mean_per_point_error.mat');
resultStdPPointSavePath = strcat(rootPath,'/result/std_per_point_error.mat');
resultCDFSavePath = strcat(rootPath,'/result/cdf.mat');
resultConfusedSavePath = strcat(rootPath,'/result/confusion.mat');

%% Find the one with hight probability
[prob, prePosition] = max(predict_label,[], 2);
test_label = test_label + 1;
confusedPosition = zeros(max(prePosition), max(test_label));

for gg = 1 : max(test_label)
    ind = find(test_label==gg);
    errorPos = prePosition(ind);
    recoo = zeros(max(prePosition), 1);
    for vv = 1 : max(prePosition)
        dd = find(errorPos == vv);
        recoo(vv) = length(dd);
    end
    confusedPosition(:,gg) = recoo;
end


%% Calculate error
totalNum = numel(test_label);
errorIPS = zeros(totalNum,1); 
for n = 1: totalNum
    cor = test_label(n);
    uncert =  prePosition(n);
    errorIPS(n) = err_table{cor}(uncert);
end

meanError = mean(errorIPS);
stdError = std(errorIPS);
resultPriMean = sprintf('The mean error of the IPS is %.4f', meanError);
resultPriStd = sprintf('The standard deviation of the IPS is %.4f', stdError);

disp(resultPriMean)
disp(resultPriStd)

%% Calculate error for each position
errPerPos = cell(nClass,1); 
meanPerPos = zeros(nClass,1);
stdPerPos = zeros(nClass,1);

for n = 1: nClass
    cor = (test_label==n);
    errPerPos{n} = errorIPS(cor);
    meanPerPos(n) = mean(errorIPS(cor));
    stdPerPos(n) = std(errorIPS(cor));
end

figure(1)
x = 1:length(meanPerPos);
stem(x, meanPerPos)
set(gca,'XTick',1:1:nClass)
ylabel('mean') 
xlabel('Position #') 

figure(2)
stem(x,stdPerPos)
set(gca,'XTick',1:1:nClass)
ylabel('std') 
xlabel('Position #') 

%% CDF
diserr = zeros(12,1);
for mm = 1:length(diserr)
    sumer = sum(errorIPS <=mm);
    diserr(mm) = sumer / totalNum;
end
diserr = [0; diserr];

x = 0:12;
figure(3)
plot(x,diserr)
set(gca,'XTick',0:1:12)
ylabel('CDF') 
xlabel('Distance error (m)') 

%% Plot confused positions
x = 1:max(test_label);
[maxs,indx] = max(confusedPosition);
% figure(4)
% pie(x, indx)
% set(gca,'XTick',1:1:max(test_label))
% ylabel('CDF') 
% xlabel('Distance error (m)') 

%% Save Results
save(resultMeanSavePath, 'meanError')
save(resultStdSavePath, 'stdError')
save(resultMeanPPointSavePath, 'meanPerPos')
save(resultStdPPointSavePath, 'stdPerPos')
save(resultCDFSavePath, 'diserr')
save(resultConfusedSavePath, 'confusedPosition')

