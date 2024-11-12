function [persist_image] = persistence(sc_image, persist_image, width, length, persist)
% Persistence (IIR Filter)
%    A low pass filter over time can be done known as Persistence.
%    Persistence is averaging the previous output image with the newly
%    calculated image.  By changing the weight of the old vs new you
%    can change the amount of persistence (or how long it takes old
%    stuff to fade in the image over time).  This filter is great for
%    remove speckle noise for slowly moving (or none moving) targets.
%    However, it causes blurring of fast moving targets.
%
%
% by ***Author***           
%
%   INPUTS
%       sc_image:     Echo Image after contrast enhancement
%       persist_image: Previous output image
%       width:        width of image
%       length:       length of image
%       persist:      amount of persistence, ranges from 0.0 for no persistence to 1.0 for maximum persistence
%
%   OUTPUTS
%       persist_image:      Echo image after persistence
%
  
