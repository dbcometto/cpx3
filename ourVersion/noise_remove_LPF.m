function [data] = noise_remove_LPF(data)
% Filters the data to remove noise and bias
%
% by Geoffrey Stentiford
%
%   INPUTS
%       data:           4 channel raw data from the phased array
%   OUTPUTS
%       data:           modified data
%

fils = load("denoise_fils.mat","-mat");

data = data - mean(data);
data = filter(fils.Num1, 1, data);
data = filter(fils.Num2, 1, data);