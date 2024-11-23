function [data_upsampled] = upsampling(data, lpf_filter_num, data_upsampled, N)
%    Upsampling by 2 is done to increase the number of beams 2x after
%    simple delay-add beamforming.
%    
%    Here, upsampling is done by 1) filling in a zeros matrix with the data
%    so that the data is zero padded, then 2) interpolating with an order-3
%    low pass filter and scaling by 2.
%
%    Note that because we only are interested in data below 10 KHz, our LPF
%    transition band is very large, which enables a very low order, and
%    thus fast, filter.
%
% by Ben Cometto
%     11/22/2024
%
%   INPUTS
%       data:           4 channel raw data from the phased array
%       lpf_filter_num: numerator coefficients for interpolating LPF
%       data_upsampled: a matrix of zeros with twice the height of data
%       N:              the height of data_upsampled
%   OUTPUTS
%       data_upsampled: upsampled data

% First, pad the data with zeros by filling in prearranged zero vector
data_upsampled(1:2:N,:) = data;

% Second, LPF to interpolate
data_upsampled = 2*filter(lpf_filter_num,1,data_upsampled);












