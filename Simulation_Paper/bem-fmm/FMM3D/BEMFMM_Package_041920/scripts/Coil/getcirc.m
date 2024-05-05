

function [x,y] = getcirc(cx, cy, radius)

theta = linspace(0,2*pi,100); %you can increase this if this isn't enough yet
x=cx + radius*cos(theta);
y=cy + radius*sin(theta);

