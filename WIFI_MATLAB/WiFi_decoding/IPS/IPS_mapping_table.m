%%
%
%  PREFORMATTED
%  TEXT
%
clear
clc

%% Map Creation
xx = 0.45 : 0.90 : 0.90 + 0.45;
yy = 0.45 : 0.90 : 0.90 * 3 + 0.45;
mapTestArray = zeros(4, 2);
for y = 1 : length(yy)
    for x = 1:length(xx)
        mapTestArray(y,x) = xx(x) + 1i*yy(y);
    end
end

xxx = 0 : 0.90 :0.90 * 2;
yyy = 0 : 0.90 : 0.90 * 4;
mapArray = zeros(5, 3);
for y = 1 : length(yyy)
    for x = 1:length(xxx)
        mapArray(y,x) = xxx(x) + 1i*yyy(y);
    end
end

%% Calculate the distance
dis_table = cell(size(mapTestArray));

for m = 1 : 4
    for k = 1: 2
        intt = mapTestArray(m, k) - mapArray;
        dis_table{m, k} = abs(intt);
    end
end

%% Actual used part in the table ( lab points)
error_table = cell(8, 1);

%% Delete the positions where test data is not covered (manul operation of choosing point)
uu = uu + 1;
error_table{uu} = dis_table{1,2};

%% Delete the positions where train data is not covered (manul operation of choosing point)
err_table = cell(43,1);
uu = uu + 1;

for vv = 1:43
    err_table{vv}(uu) = error_table_2{vv}(1,1);
end

%%
err_table = cell(8,1);
for i = 1 : 8
     kkk = error_table{i,:} ;
     ttt = flipud(kkk);
     ccc = reshape(ttt, 15, 1);
     err_table{i} = ccc;
end