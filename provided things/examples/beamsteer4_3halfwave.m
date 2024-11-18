% simulates a 4 channel array spaced at 1.5 wavelengths
% Note the ambiguity due to multiple constructive inteference
% (or multiple beams)
lambda = 1;
r0 = 5;
d = 3*lambda/2;
xoff = d/2;




for r = 0:.1:r0
    x = -5:.001:5;
    y1 = sqrt(r.^2-(x+xoff).^2);
    y2 = sqrt(r.^2-(x-xoff).^2); 
    y3 =  sqrt((r-lambda).^2-(x+xoff).^2);
    y4 =  sqrt((r-lambda).^2-(x-xoff).^2);
    h = plot(x,y1,x,y2,x,y3,x,y4);
    axis([-5 5 0 10])
    drawnow;
end