function [data] = cal_start(data, N, pulse_length)
% Detects start of transmitted pulse, when signal rises over "threshold",
% then blanks out the transmitted pulse, shifts the data to align
% "time zero" at the start of the array, and blanks out the
% post-reflection noise
%
% by Geoffrey Stentiford
%
%   INPUTS
%       data:           4 channel raw data from the phased array
%       N:              The number of elements in each array
%       pulse_length:   the number of samples in the transmitted pulse
%   OUTPUTS
%       data:           modified data
%       data_start:     location in data at end of blanking region
%       data_end:       location in data at the end of "good" data

temp = abs(data);
count = 0;
idx = 1;

while count < 10 && count < length(data(:,1))
    if temp(idx,1) > .08
        count = count + 1;
    else
        count = 0;
    end
    idx = idx + 1;
end
idx = idx - 1;
begin = idx;
count = 0;
while count < 30 && count < length(data(:,1))
    if temp(idx,1) < .04
        count = count + 1;
    else
        count = 0;
    end
    idx = idx + 1;
end
idx = idx - 1;

data(1:1:idx,:) = data(1:1);
%data(fin:1:end,:) = data(1:1);
data(1:1:end-begin,:) = data(begin+1:1:end,:);
data(1+end-begin:1:end,:) = data(1:1);

his = max(temp(idx:1:end,:));
rels = his/max(his);

for n = 1:1:4
    data(:,n) = data(:,n)/rels(n);
end