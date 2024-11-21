function [data_upsampled] = upsampling(data, upsample_by)
%    Upsampling is done to increase the number of beams after
%    beamforming (assuming simple delay-add beamforming)
%    Assuming signal is bandlimited (by stage 5), upsampling
%    can be done my (1) padding with zeroes between samples
%    (2) then low pass filter to remove the signal "images" created
%
% by ***Author***
%
%   INPUTS
%       data:           4 channel raw data from the phased array
%       upsample_by:    upsample rate (currently only supports "2")
%   OUTPUTS
%       data_upsampled:           modified data
%








