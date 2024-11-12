% simulates a 4 channel array spaced at 1/2 wavelength
% Note the constructive inteference only occurs along one path
% (or beam) at a time

r0 = 10;
x0 = 0.5;
d = 10;
xmin = -d;
xmax = d;
ymin = 0;
ymax = d;
max_const = d/3;  %set this to where the max construction will occur
for theta = -90:10:90
x = xmin:.01:xmax;
    ro1 = sqrt((max_const)^2 + (3*x0)^2 + 6*x0*max_const*cos((90-theta)*pi/180))-max_const;
    ro2 = sqrt((max_const)^2 + (x0)^2 + 2*x0*max_const*cos((90-theta)*pi/180))-max_const;
    ro3 = sqrt((max_const)^2 + (x0)^2 - 2*x0*max_const*cos((90-theta)*pi/180))-max_const;
    ro4 = sqrt((max_const)^2 + (3*x0)^2 - 6*x0*max_const*cos((90-theta)*pi/180))-max_const;
for r = 0:.1:r0
    y1 = sqrt((r+ro1).^2 - (x+3*x0).^2);
    y2 = sqrt((r+ro2).^2 - (x+x0).^2);
    y3 = sqrt((r+ro3).^2 - (x-x0).^2);
    y4 = sqrt((r+ro4).^2 - (x-3*x0).^2);
    r1 = 0:0.1:r;

    xsol = .5*((ro2)-(ro3)).*(r1 + 0.5*((ro2)+(ro3)))/x0;
    ysol = sqrt((r1.^2)-(xsol).^2);
    h = plot(x,y1,x,y2,x,y3,x,y4,xsol,ysol);
    axis([xmin xmax ymin ymax])
    drawnow;
end
end
