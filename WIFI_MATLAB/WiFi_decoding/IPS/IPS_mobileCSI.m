%%
% 
%  Observe dynamic csi data
% 
% 
% clear
% close all
% clc
%%
% csiData6 = reshape(dataInput.', numel(dataInput), 1);
% testData= [csiData1,csiData2,csiData3,csiData4,csiData5,csiData6];
trainData= [csiData1,csiData2,csiData3,csiData4,csiData5,csiData6];
%%
yyy = 4000;
xxx = 1000000;
for i = 1 : 6
    subplot(6, 1, i)
    plot(trainData(1:xxx,i))
    xlabel("sample")
    ylabel("amplitude")
%     title(['point ', num2str(i+4)])
    ylim([0 yyy])
end
%%
x1 = trainData(1:510000,1);
x2 = trainData(16000:530000,2);
x3 = trainData(1:92500,3);
x3 = [x3;trainData(193400:end,3)];
x4 = trainData(1:383200,4);
x4 =[x4; trainData(475800:end,4)];
x5 = trainData(1:510000,5);
x6 = trainData(1:242000,6);
x6 = [x6;trainData(320000:end,6)];
%%
% x_ = [x1,x2,x3,x4,x5,x6];
    subplot(6, 1, 1)
    plot(x1)
%     title(['point ', num2str(i+4)])
   xlim([0 510000])
    ylim([0 yyy])
        subplot(6, 1, 2)
    plot(x2)
%     title(['point ', num2str(i+4)])
   xlim([0 510000])
    ylim([0 yyy])
        subplot(6, 1, 3)
    plot(x3)
%     title(['point ', num2str(i+4)])
   xlim([0 510000])
    ylim([0 yyy])
        subplot(6, 1, 4)
    plot(x4)
%     title(['point ', num2str(i+4)])
   xlim([0 510000])
    ylim([0 yyy])
        subplot(6, 1, 5)
    plot(x5)
%     title(['point ', num2str(i+4)])
   xlim([0 510000])
    ylim([0 yyy])
        subplot(6, 1, 6)
    plot(x6)
%     title(['point ', num2str(i+4)])
   xlim([0 510000])
    ylim([0 yyy])
    


%%
close all
csiData1_reshape = reshape(csiData1,56,numel(csiData1)/56);
figure(1)
for i = 1 : 56
    % dbscan
    [denoiseCsi, epss, numOutlier] = denoise_dbscan(csiData1_reshape(i, 2000:3400), 50);
%     plot(csiData1_reshape(i,2000:3400))
    plot(denoiseCsi)
    hold on
        ylim([0 3500])
end
csiData2_reshape = reshape(csiData2,56,numel(csiData1)/56);
figure(2)
for i = 1 : 56
    plot(csiData2_reshape(i,4500:5900))
    hold on
        ylim([0 3500])
end
csiData3_reshape = reshape(csiData3,56,numel(csiData1)/56);
figure(3)
for i = 1 : 56
    plot(csiData3_reshape(i,5500:7400))
    hold on
        ylim([0 3500])
end

