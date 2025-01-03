function [sc_image] = contrast(sc_image, bottom, top)
% Constrast/Brightness
%   This function warp's the dynamic range of the image intensity
%   to emphasis features of interest and minimize noise.
%   There are many ways to do this such as:
%          - logrithmic compression:  I = exp(I)
%          - "window/level" (known as imadjust in matlab)
%          - histeq (in matlab)
%  
%
% by Geoffrey Stentiford
%
%   INPUTS
%       sc_image:     Echo Image after scan conversion
%       width:        width of image
%       length:       length of image
%       top:          top threshold for window/level (must be between 0 and 1) 
%       bottom:       bottom thrseshold for window/level (must be between 0 and 1)
%
%   OUTPUTS
%       sc_image:      Echo image after contrast enhancement
%
  

high = max(max(sc_image));
sc_image = contrastFix(sc_image, top, bottom, high);