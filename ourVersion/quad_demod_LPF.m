function [demod_I_LPF, demod_Q_LPF] = quad_demod_LPF(demod_I, demod_Q, NumBeams, filter_coef, demod_I_LPF, demod_Q_LPF)
%  LPF to remove harmonic at twice f, keeping DC component
%  Convolves data with Filter_coef of length WindowLength
%
% by Victor Chen and Geoffrey Stentiford
%
%   INPUTS
%       demod_I:        Inphase component of the mixed beamformed data
%       demod_Q:        quadrature component of the mixed beamformed data
%       NumBeams:       number of beams
%       FrameSize:      Number of samples in each beam
%       WindowLength:   WindowLength of filter
%       filter_coef:    filter coefficients to convolve with
%
%   OUTPUTS
%       demod_I:        Inphase component after LPF
%       demod_Q:        quadrature component after LPF
%

% TEST DATA - UNCOMMENT TO USE
% Using a sinwave for all channels, we can test the beamform on various
% frequencies to see waves in our final beamform image.
% x = 2    % adjust x for different frequencies
%t = 1:4000
%channel = sin(t/x)
%data2 = [channel', channel', channel', channel']

% Loop through each signal
for i = 1:NumBeams
    % Convolve the current signal with the filter
    demod_I_LPF(i, :) = conv(demod_I(i, :), filter_coef, 'same');
    demod_Q_LPF(i, :) = conv(demod_Q(i, :), filter_coef, 'same');
end


