function [demod_I_LPF, demod_Q_LPF] = quad_demod_LPF(demod_I, demod_Q, NumBeams, FrameSize, WindowLength, filter_coef)
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
%       demod_I_LPF:        Inphase component after LPF
%       demod_Q_LPF:        quadrature component after LPF
%



% Apply the LPF
demod_I = filter(filter_coef, 1, demod_I);
demod_Q = filter(filter_coef, 1, demod_Q);

demod_I_LPF = demod_I
demod_Q_LPF = demod_Q


