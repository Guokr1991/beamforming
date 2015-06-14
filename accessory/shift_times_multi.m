function [pad_rf, s0, t0] = shift_times_multi(rf, st, fs)
% Written by WJL - 10/6/14
% Parameters:
%     rf - cell array of matrices of rf at each receive channel for all tx
%     events
%     st - array containing the shifted start times of each transmit event
%     fs - sampling frequency used for acquisition
% Returns: 
%     pad_rf - array of padded rf matrices with same t0, 3rd dim 
%     indicates each tx
%     t0 - start time of padded rf (given in sample number)

for j = 1:length(st)
    tot_sample(j) = st(j)*fs+size(rf{j},1);
end
min_sample = min(st)*fs;
max_sample = max(tot_sample);

if min_sample ~= round(min_sample)
    disp(['Calculated min sample: ' num2str(min_sample)]);
    disp(['Rounded min sample: ' num2str(round(min_sample))]);
    min_sample = round(min(st)*fs);

elseif max_sample ~= round(max_sample)
    disp(['Calculated max sample: ' num2str(max_sample)]);
    disp(['Rounded max sample: ' num2str(round(max_sample))]);
    max_sample = round(max_sample);
end

n_rcv = size(rf{1},2);
n_tx = length(st); 
pad_rf = zeros(max_sample-min_sample, n_rcv, n_tx);
for j = 1:n_tx
    tmp_rf = rf{j};
    if st(j)*fs ~= min_sample
        tmp_rf = [zeros(round(st(j)*fs) - min_sample, n_rcv); tmp_rf];
    end
    if tot_sample(j) ~= max_sample
        tmp_rf = [tmp_rf; zeros(max_sample - tot_sample(j), n_rcv)];

    end
    pad_rf(:,:,j) = tmp_rf;
    clear tmp_rf
end

s0 = min_sample;
t0 = min(st);


