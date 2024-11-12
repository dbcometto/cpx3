function [beams] = beamform(data2, NumBeams, FrameSize, num_elements)
%    Convert the data from the 4 channels into 21 "beams" using simple
%    delay-sum beamforming
%
%
% by ***AUTHOR****
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





