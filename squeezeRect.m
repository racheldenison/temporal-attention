function newRect = squeezeRect(rect, squeezeXY)
%
% function newRect = squeezeRect(rect, squeezeXY)
%
% rect is a Psychtoolbox rect, a vector with 4 elements
% Note a PTB rect is of the form [left top right bottom]
%
% squeezeXY is the proportion of the original dimension by which to squeeze
% a dimension. 1 or 2 element vector [squeezeX squeezeY]. If just 1 
% element, will use the same value for both dimensions. Values less than
% 1 shrink the dimension, values greater than 1 expand the dimension.
%
% Rachel Denison
% March 2014

% deal with inputs
if length(squeezeXY)==1
    squeezeXY = [squeezeXY squeezeXY];
elseif length(squeezeXY)>2
    error('squeezeXY must have exactly 1 or 2 elements')
end

% squeeze
xSz = rect(3)-rect(1);
ySz = rect(4)-rect(2);

cx = mean(rect([1 3]));
cy = mean(rect([2 4]));

newXSz = xSz*squeezeXY(1);
newYSz = ySz*squeezeXY(2);

newRect = CenterRectOnPoint([0 0 newXSz newYSz], cx, cy);