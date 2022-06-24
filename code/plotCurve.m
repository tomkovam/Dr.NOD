function plotCurve(x1, x2, y1, y2, p, colour, lineWidth)
% p = 1: plots half a circle
% p < 1: plots the top p-fraction of the half a circle

if (x2 < x1)
    tmp = x2;
    x2 = x1;
    x1 = tmp;
end
if (y2 < y1)
    tmp = y2;
    y2 = y1;
    y1 = tmp;
end


% figure(10);
% % th = linspace(pi/2, -pi/2, 100);
% th = linspace(-pi/2, pi/2, 100);
% R = 1;  %or whatever radius you want
% x = R*cos(th) + 5;
% y = R*sin(th) + 4;
% plot(x,y);

% function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)

rX = (x2 - x1)/2;
rY = (y2 - y1)/2;
xMiddle = x1 + rX;
yBottom = y1;
n = 1000;
k = max([1, round(n*(1-p))]);


th = linspace(0, pi, n);
ang = th(k:end-k+1);
xp=rX*cos(ang);
yp=rY*sin(ang);
yp = yp - min(yp);
plot(xMiddle+xp,yBottom+yp, 'Color', colour, 'LineWidth', lineWidth);
