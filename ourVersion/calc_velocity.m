function velocity = calc_velocity([row1,col1], [row2,col2], fpsRow, fpsCol, Fs)

% by C1C Nicholas Csicsila
% 
%   Uses 2-D correlation function to correlate with each beam before
%   demodulation, in order to find 
%   INPUTS
%       [row1, col1]:   location of first track
%       [row2, col2]:   location of second track
%       fpsRow:         Feet per second for Row
%       fpsCol:         Feet per second for Column
%       Fs:             Sampling frequency    
%   OUTPUTS
%       velocity:       returns single velocity value in ft/s

time = 1/Fs;

distance = sqrt(((col2 - col1)*fpsCol)^2 + ((row2 - row1)*fpsRow)^2))

velocity = distance/time