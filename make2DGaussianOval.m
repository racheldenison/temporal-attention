function gaussian = make2DGaussianOval(w, h, x0, y0, ... 
    sigmaX, sigmaY, gaussAmp)
%
% function gaussian = make2DGaussianOval(w, h, x0, y0, ...
%   sigmaX, sigmaY, gaussAmp)


[x y] = meshgrid(1:w, 1:h);

gaussian = gaussAmp*exp(-(((x-x0)/(2*sigmaX)).^2)-(((y-y0)/(2*sigmaY)).^2));
