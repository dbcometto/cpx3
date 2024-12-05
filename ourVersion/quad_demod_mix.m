function [demod_I, demod_Q] = quad_demod_mix(beams, cos_table, sin_table)
% Demodulate the echo image by finding the envelope of the carrier pulse
% First, multiply by the carrier (cos + jsin)
% This moves the echo signal to DC and makes another copy at twice the
% carrier frequency.
%
% The cos and sin table are assumed precomputed to speed up the computation
%

% by Geoffrey Stentiford
%
%   INPUTS
%       beams:          beamformed data
%       cos_table:      pre-computed lookup table for cosine
%       sin_table:      pre-computed lookup table for sine
%
%   OUTPUTS
%       demod_I:        Inphase component of the mixed beamformed data
%       demod_Q:        quadrature component of the mixed beamformed data
%

dims = size(beams);
demod_I = zeros(dims(1), dims(2));
demod_Q = zeros(dims(1), dims(2));
for i = 1:dims(2)
    demod_I(:,i) = beams(:,i) .* cos_table(i);
    demod_Q(:,i) = beams(:,i) .* sin_table(i);
    %demod_I(:,i) = beams(:,i) * cos(factor*i);
    %demod_Q(:,i) = beams(:,i) * sin(factor*i);
end
