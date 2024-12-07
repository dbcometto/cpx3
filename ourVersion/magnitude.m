function [Mag_image] = magnitude(demod_I_LPF, demod_Q_LPF)
%  Find magnitude (echo image) from imaginary
%           Find magnitude of I + jQ
%  
% by Victor Chen and Geoffrey Stentiford
%
%   INPUTS
%       demod_I_LPF:    Inphase component of the mixed/LPF beamformed data
%       demod_Q_LPF:    quadrature component of the mixed/LPF beamformed data
%       NumBeams:       number of beams
%       FrameSize:      Number of samples in each beam
%
%   OUTPUTS
%       Mag_image:      Echo image after magnitude
%

Mag_image = sqrt(demod_I_LPF.^2 + demod_Q_LPF.^2);