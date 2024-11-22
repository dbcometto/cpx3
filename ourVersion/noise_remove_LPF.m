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

load denoise_fils.mat -mat

parfor i = 1:4
    data(:,i) = data(:,i) - mean(data(:,i));
    data(:,i) = filter(Num1, 1, data(:,i));
    data(:,i) = filter(Num2, 1, data(:,i));
end