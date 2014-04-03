function t = makeUprightTriangle(imW, imH, x0, y0, triW, triH)
%
% t = makeUprightTriangle(imW, imH, x0, y0, triW, triH)
%
% INPUTS:
% imW and imH are the width and height of the image
% x0 and y0 are the coordinates of the center of the triangle
% triW and triH are the width and height of the triangle
%
% OUTPUTS:
% t is an [imW x imH] image with a [triW x triH] triangle centered at (x0,
% y0). the area in the triangle = 1 and the background = 0.

if nargin==0
    imH = 200;
    imW = 300;
    x0 = 200;
    y0 = 50;
    triW = 50;
    triH = 100;
end

rect = CenterRectOnPoint([0 0 triW triH], x0, y0); % [left top right bottom]

x = [rect(1) mean(rect([1 3])) rect(3)]; % [bottom left, top, bottom right] 
y = [rect(4) rect(2) rect(4)];
t = poly2mask(x, y, imH, imW);

t = double(t);