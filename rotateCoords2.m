function coordsRot = rotateCoords2(coords, rotation)
%
% function coordsRot = rotateCoords2(coords, rotation)
%
% coords have 2 rows, x and y
% rotation angle in degrees, positive is clockwise
%
% unlike rotateCoords, does not assume that the coords are centered at
% (0,0). will center the object based on the mean of the x and y coords for
% the purposes of the rotation, and then shift the object back to its
% original position.
%
% Modified from rotateCoords.m
%
% Rachel Denison
% March 2015

rotation = rotation*pi/180;

x = coords(1,:);
y = coords(2,:);

dcX = mean(x);
dcY = mean(y);

x = x-dcX;
y = y-dcY;

r = sqrt(x.^2+y.^2);
theta = atan2(y,x);

% apply rotation
% unit circle goes counterclockwise, but we want clockwise rotation, so
% subtract
thetaRot = theta - rotation; 

xRot = r.*cos(thetaRot);
yRot = r.*sin(thetaRot);
% [y,x] = pol2cart(theta,r)

coordsRot(1,:) = xRot + dcX;
coordsRot(2,:) = -yRot + dcY; % y is negative in PTB
