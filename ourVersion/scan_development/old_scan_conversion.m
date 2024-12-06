function [image] = scan_conversion(Mag_image, ind_bkn, ind_bk1n, ind_bkn1, ind_bk1n1, BMAM, BMA, BAM, BA)
% Scan Conversion
%   Uses complex lookup tables to properly plot the "beams" in their 
%   actual spatial geometry.
%   Implements polar conversion from rectangular array, and bilinear
%   interpolation of  pixels that lie in between beam lines.
%
% by Ben Cometto
%
%   INPUTS
%       Mag_image   :     image before scan conversion
%       ind_bkn     :     indices for beam data at k and n
%       ind_bk1n    :     indices for beam data at k+1 and n
%       ind_bkn1    :     indices for beam data at k and n+1
%       ind_bk1n1   :     indices for beam data at k+1 and n+1
%       BMAM        :     per pixel (1-beta)*(1-alpha) values [fractional distance between beams]
%       BMA         :     per pixel (1-beta)*(alpha) values
%       BAM         :     per pixel (beta)*(1-alpha) values
%       BA          :     per pixel (beta)*(alpha) values
%
%   OUTPUTS
%       image       :     scan converted image
%
  
image=BMAM.*Mag_image(ind_bkn)+BMA.*Mag_image(ind_bk1n) + BAM.*Mag_image(ind_bkn1) + BA.*Mag_image(ind_bk1n1);