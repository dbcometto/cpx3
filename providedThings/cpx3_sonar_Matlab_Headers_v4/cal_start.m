function [data, data_start, data_end] = cal_start(data, N, pulse_length)
% Detects start of transmitted pulse, when signal rises over "threshold",
% then blanks out the transmitted pulse, and shifts the data to aligne
% "time zero" at the start of the array.
%
% by ***Author***
%
%   INPUTS
%       data:           4 channel raw data from the phased array
%       N:              The number of elements in each array
%       pulse_length:   the number of samples in the transmitted pulse
%   OUTPUTS
%       data:           modified data
%       data_start:     location in data at end of blanking region
%       data_end:       location in data at the end of "good" data