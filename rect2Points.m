function points = rect2Points(rect)
%
% function points = rect2Points(rect)
%
% converts a PTB rect to a list of xy points (2 x n points)
% a rect is [upperleftX upperleftY lowerrightX lowerrightY]

% upper left
points(1,1) = rect(1);
points(2,1) = rect(2);

% upper right
points(1,2) = rect(3);
points(2,2) = rect(2);

% lower right
points(1,3) = rect(3);
points(2,3) = rect(4);

% lower left
points(1,4) = rect(1);
points(2,4) = rect(4);