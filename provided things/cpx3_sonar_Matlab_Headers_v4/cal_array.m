function [data] = cal_array(data, data_start, data_end)
% Self-Calibrates a 4 channel array
% - removes DC bias
% - normalizes energy across arrays
%
% by ***Author***
%
%   INPUTS
%       data:           4 channel raw data from the phased array
%       data_start:     location in data at end of blanking region
%       data_end:       location in data at the end of "good" data
%   OUTPUTS
%       data:           modified data



