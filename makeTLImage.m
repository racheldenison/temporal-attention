function im = makeTLImage(sz, lineWidth, propOffset, color, bgColor)

% function im = makeTLImage(sz, lineWidth, propOffset, color, bgColor)
%
% Makes a T/L (eg. visual search) image
% When propOffset = 0, it is a T. When propOffset > 0, the trunk of the T
% is shifted to the right. When propOffset < 0, the trunk of the T is
% shifted to the left.
%
% INPUTS:
% sz is the size of one side of the square, in pixels
% lineWidth is the width of the arms of the T/L, in pixels
% propOffset is how much the trunk of the T is offset from the
% center, as a proportion of the distance between the center and the edge
% (sz/2). if >1 or <-1, the trunk will just display on the edge of the
% image.
% color is the grayscale value of the T/L
% bgColor is the grayscale value of the background

% example inputs
% sz = 100;
% lineWidth = 5;
% propOffset = .5;
% color = 1;
% bgColor = 0.5;

% base image in background color
im = ones(sz)*bgColor;

% top of the T
im(1:lineWidth,:) = color;

% trunk of the T
xStart = sz/2 - lineWidth/2 + sz/2*propOffset;
trunkX = (1:lineWidth) + round(xStart);

% if some of the trunk is out of the picture, just move it right back in.
% this is a bit of a hack, but it means the two arms of the T/L will always
% be the same width
spillover = trunkX(end) - sz;
if spillover>0
    trunkX = trunkX - spillover;
end
if trunkX(1)<1
    trunkX = trunkX + 1-trunkX(1);
end

% fill in the trunk
im(:,trunkX) = color;

% show the image
% imshow(im)