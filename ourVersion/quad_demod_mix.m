function [demod_I, demod_Q] = quad_demod_mix(beams, NumBeams, FrameSize, frequency, SampleRate)
% Demodulate the echo image by finding the envelope of the carrier pulse
% First, multiply by the carrier (cos + jsin)
% This moves the echo signal to DC and makes another copy at twice the
% carrier frequency.
%
% The cos and sin table are assumed precomputed to speed up the computation
%

% by Victor Chen and Geoffrey Stentiford
%
%   INPUTS
%       beams:          beamformed data
%       NumBeams:       number of beams
%       FrameSize:      Number of samples in each beam
%       Do you want to add some precomputed look-up-tables?
%
%   OUTPUTS
%       demod_I:        Inphase component of the mixed beamformed data
%       demod_Q:        quadrature component of the mixed beamformed data
%

factor = pi*.1;
dims = size(beams);
demod_I = zeros(dims(1), dims(2));
demod_Q = zeros(dims(1), dims(2));
for i = 1:dims(2)
    demod_I(:,i) = beams(:,i) * cos(factor*i);
    demod_Q(:,i) = beams(:,i) * sin(factor*i);
end
