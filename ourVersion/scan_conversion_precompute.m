function [ind_bkn, ind_bk1n, ind_bkn1, ind_bk1n1, BMAM, BMA, BAM, BA, rppf, cppf] = scan_conversion_precompute(fc, Fs, C, h, w, shift, imager, imagec)
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


ishiftc = floor(imagec/2);


% Precompute thetas
thetas = zeros(h,w);
for kk=1:w
    k = kk-shift;
    thetas(:,kk) = asin(2*k*fc/Fs);
end

% Precompute rs
rs = zeros(h,w);
for n=0:h-1
    rs(n+1,:) = n*C/Fs/2;
end
rs;

% Precompute max distance and ratios
max_distance = 4000*C/Fs/2;

rppf = (imager-1)/rs(end,1); % pixels per foot
cppf = (imagec/2-1)/rs(end,1);

% Precomute indexes and alphas/betas
KKNvals = zeros(imager,imagec,2);
KK1N1vals = zeros(imager,imagec,2);
Alphas = zeros(imager,imagec);
Betas = zeros(imager,imagec);

for pc = 1:imagec
    for pr = 1:imager

        % Establish x and y
        px = (pc-1-ishiftc)/cppf;
        py = (pr-1)/rppf;        

        % Establish rp and thetap
        rp = sqrt(px^2+py^2);
        thetap = atan(px/py);

        % Find indexes (kk and n)
        if (px ~= 0) % handle edge case
            kp = sin(thetap)/2*Fs/fc;
        else
            kp =  0;
        end
          
        k = floor(kp);
        kk = k+shift;

        np = rp/C*Fs*2;
        n = floor(np)+1;  

        if (kk < 21) % bound kk
            KKNvals(pr,pc,1) = kk;
            KK1N1vals(pr,pc,1) = kk+1;
        elseif (kk==21)
            KKNvals(pr,pc,1) = 21;
            KK1N1vals(pr,pc,1) = 21;
        else
            KKNvals(pr,pc,1) = 1;
            KKNvals(pr,pc,2) = 1;
            KK1N1vals(pr,pc,2) = 1;
            KK1N1vals(pr,pc,2) = 1;
        end

        if (n < 4000) % bound n
            KKNvals(pr,pc,2) = n;
            KK1N1vals(pr,pc,2) = n+1;
        elseif (n==4000)
            KKNvals(pr,pc,2) = 4000;
            KK1N1vals(pr,pc,2) = 4000;
        else
            KKNvals(pr,pc,1) = 1;
            KKNvals(pr,pc,2) = 1;
            KK1N1vals(pr,pc,1) = 1;
            KK1N1vals(pr,pc,2) = 1;
        end

        % Then, calculate alphas and betas 
        if (n>=4000)||(kk>=21)||(kp==k)||(np==n)
            alpha = 0;
            beta = 0;
        else
            r = rs(n,1);
            r1 = rs(n+1,1);
            dr = r1-r;
    
            theta = thetas(1,kk);
            theta1 = thetas(1,kk+1);
            dtheta = theta1-theta;
            
            alpha = (thetap - theta)/dtheta;
            beta = (rp - r)/dr;  
        end

        Alphas(pr,pc) = alpha;
        Betas(pr,pc) = beta;

    end
end


% Precompute linear indexes
ind_bkn = sub2ind([21 4000],KKNvals(:,:,1),KKNvals(:,:,2));
ind_bk1n = sub2ind([21 4000],KK1N1vals(:,:,1),KKNvals(:,:,2));
ind_bkn1 = sub2ind([21 4000],KKNvals(:,:,1),KK1N1vals(:,:,2));
ind_bk1n1 = sub2ind([21 4000],KK1N1vals(:,:,1),KK1N1vals(:,:,2));


% Precompute Scaling Factors (Beta/Alpha, Minus)
BMAM = (1-Betas).*(1-Alphas);
BMA = (1-Betas).*Alphas;
BAM = Betas.*(1-Alphas);
BA = Betas.*Alphas;