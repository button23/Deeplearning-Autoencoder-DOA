function [dataOutput, epss, numOutlier] = denoise_dbscan(dataInput, inEps)
numOutlier = zeros(1000, 1);
epss = zeros(1000, 1);
for i = 1 : 1
    close all;
    % the parameter for dbscan
    eps = inEps;
    % apply dbscan algorithm
    chEst_lltf_ind = dbscan(dataInput, eps, 30);
    discardd = sum(chEst_lltf_ind ~= 1);
    while true
        if discardd >= numDiscard
            eps = eps + 1;
            chEst_lltf_ind = dbscan(dataInput, eps, 30);
            discardd = sum(chEst_lltf_ind ~= 1);
        else
            break;
        end
    end
    % recorder the number of discarded samples
    epss(i) = eps;
    numOutlier(i) = discardd;
    
%     figure(i)
%     hold on
    indd = find(chEst_lltf_ind == 1);
    dataOutput = dataInput(indd, :);
%     for k = 1: length(indd)
%         plot(dataOutput(k, :));
%     end
end
end