function [data] = time_gain_compensation(data, timeGainValue)

% Adjusts for annenuation of signal over distance
%  linear or log?     tested with tone at proper distance
%  For now: APPEARS to be linear; slope = falls 1.0 over 600 samples
%
% by C1C Nicholas Csicsila
%
%   INPUTS
%       data:           4 channel raw data from the phased array
%       N:              number of samples in each channel
%   OUTPUTS
%       data:           modified data
   
data = (data .* timeGainValue);











