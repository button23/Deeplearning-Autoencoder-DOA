%%
num_loc = 5;
doa = zeros(4,num_loc);
for i = 1 :4 % 4 APs
    for j = 1 :num_loc % the total number of locations
        mid = rays{i,j};

        if isempty(mid)
            doa(i,j) = 0;
        else
            doa(i,j) = mid(1).AngleOfArrival(1);
        end
    end
end

%% add 180 to make all the angles to be positive
doa_shift_180 = doa +180;
doa_shift_180_reshape = reshape(doa_shift_180,num_loc*4,1);
min(doa_shift_180_reshape)
doa_input = sum(doa_shift_180, 1);
doa_input = int16(doa_input);
doa_input = doa_input - min(doa_input);

uniquevalue = unique(doa_input);
save('doa_input.mat','doa_input')

%% add 180 to make all the angles to be positive
doa_csi_ = reshape(doa_csi, 4, num_loc);
doa_shift_180_2 = doa_csi_ +180;
doa_shift_180_reshape_2 = reshape(doa_shift_180_2,num_loc*4,1);
min(doa_shift_180_reshape_2)
doa_input_ = sum(doa_shift_180_2, 1);
doa_input_ = int16(doa_input_);
doa_input_ = doa_input_ - min(doa_input_);

uniquevalue = unique(doa_input_);
save('doa_input.mat','doa_input')