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

for i = 1:NumBeams
    
end
    




