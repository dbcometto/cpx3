function [demod_I, demod_Q] = quad_demod_mix(beams, NumBeams, cos_table, sin_table)
% Demodulate the echo image by finding the envelope of the carrier pulse
% First, multiply by the carrier (cos + jsin)
% This moves the echo signal to DC and makes another copy at twice the
% carrier frequency.
%
% The cos and sin table are assumed precomputed to speed up the computation

% by Geoffrey Stentiford
%
%   INPUTS
%       beams:          beamformed data
%       NumBeams:       number of beams
%       cos_table:      pre-computed lookup table for cosine
%       sin_table:      pre-computed lookup table for sine
%
%   OUTPUTS
%       demod_I:        Inphase component of the mixed beamformed data
%       demod_Q:        quadrature component of the mixed beamformed data
%

demod_I = beams .* cos_table(ones(NumBeams, 1), :);
demod_Q = beams .* sin_table(ones(NumBeams, 1), :);
