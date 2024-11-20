function [beams] = beamform(data2, NumBeams, FrameSize, num_elements)
%    Convert the data from the 4 channels into 21 "beams" using simple
%    delay-sum beamforming
%
%
% by VICTOR CHEN
%
%   INPUTS
%       data2:          4 channel data from the phased array
%       NumBeams:       number of beams = ratio of Fs/F + 1
%                       we upsampled to increase the number of beams
%       FrameSize:      Number of samples in each channel (after upsample)
%       num_elements:   Number of channels in the phased array
%
%   OUTPUTS
%       beams:           beamformed data
%
 
% Figure out format of data2, may need to reshape.
 
% Initialize n and k terms
n = 0;
k = 0;
 
% Create array to hold beam vals
 
% Create function for adding all channel sample data to create beam.
for i = 1:4
    X_val = data2[..] + data2[..] + ...
end
 
 
 
 


