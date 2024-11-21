function [data_upsampled] = upsampling(data, lpf_filter_num)
%    Upsampling is done to increase the number of beams after
%    beamforming (assuming simple delay-add beamforming)
%    Assuming signal is bandlimited (by stage 5), upsampling
%    can be done my (1) padding with zeroes between samples
%    (2) then low pass filter to remove the signal "images" created
%
% by Ben Cometto
%
%   INPUTS
%       data:           4 channel raw data from the phased array
%       upsample_by:    upsample rate (currently only supports "2")
%   OUTPUTS
%       data_upsampled:           modified data

% First, pad the data with zeros
data_upsampled = zeros(2*height(data),width(data));
data_upsampled(1:2:height(data_upsampled),:) = data;

% Second, LPF to interpolate
data_upsampled = 2*filter(lpf_filter_num,1,data_upsampled);












