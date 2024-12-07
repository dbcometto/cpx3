function velocity = calc_velocity(row1, col1, row2, col2, ppfRow, ppfCol, Fs, samplesApart)

% by C1C Nicholas Csicsila
% 
%   Uses 2-D correlation function to correlate with each beam before
%   demodulation, in order to find 
%   INPUTS
%       [row1, col1]:   location of first track
%       [row2, col2]:   location of second track
%       ppfRow:         Pixel per foot Row
%       ppfCol:         Pixel per foot for Column
%       Fs:             Sampling frequency
%       samplesApart:   Images between samples for distance
%   OUTPUTS
%       velocity:       returns single velocity value in ft/s

time = 50*(4000/Fs);

distance = sqrt(((col2 - col1) / ppfCol)^2 + ((row2 - row1) / ppfRow)^2);

velocity = distance/time;
