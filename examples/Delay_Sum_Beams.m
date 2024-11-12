function Delay_Sum_Beams()
% Delay_Sum_Beams.m
%    plots 11 beams (10x oversampled) for the delay-sum beamforming
%    algorithm.
% By George York    

white = 1.0;
black = 0.0;
fc = 10000;
Fs = 100000;
ratio = Fs/fc;
width = 200;
height = 100;
range = height - 2;

beam_image(1:height,1:width) = white;

for k = -(ratio/2):(ratio/2)   % for each beam
    angle = asin(2*k*(1/ratio));
    del_x = sin(angle);
    del_y = cos(angle);
    x = width/2;
    y = 1;
    for i=1:range
        beam_image(round(y+i*del_y), round(x+i*del_x)) = black;
    end
    colormap(gray)
    axis equal
    figure(1)
    imagesc(beam_image);
    pause(0.5);
end