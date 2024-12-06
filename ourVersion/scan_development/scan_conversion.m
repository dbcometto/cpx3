function [sc_image] = scan_conversion(Mag_image, NumBeams, FrameSize, min_range, max_range)
% Scan Conversion
%   Properly plot the "beams" in their actual spatial geometry
%   Requires a polar conversion (beam data is stored in a rectangular
%   array, but must be plotted polar).  Also requires interpolating
%   pixels that lie in between beam lines.
%
%   This can be greatly speeded up by precomputing lookup tables
%
%
% by ***Author***
%
%   INPUTS
%       Mag_image:    Echo Image before scan conversion
%       NumBeams:     number of beams
%       FrameSize:    Number of samples in each beam
%       max_range:    range of the output display in # of pixels
%       min_range:    range blanked on the output display in # of pixels
%
%   OUTPUTS
%       sc_image:      Scan Converted Echo image
%
  

     


