function [row, col] = tracking(data, template)
%
% by C1C Nicholas Csicsila
% 
%   Uses 2-D correlation function to correlate with each beam before
%   demodulation, in order to find 
%   INPUTS
%       data:           takes most recent array of sonar data
%       template:       can be adjusted to most accuratley fit image     
%   OUTPUTS
%       row:           row of max value
%       col:           col of max value
   
c = xcorr2(data, template);


% Find the maximum value and its linear index
[maxValue, linearIndex] = max(c(:));

% Convert linear index to row and column indices
[row, col] = ind2sub(size(c), linearIndex);
