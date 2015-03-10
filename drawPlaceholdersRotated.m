function drawPlaceholdersRotated(window, color, bgColor, rect, penWidth, on, rotation)

% function drawPlaceholdersRotated(window, color, bgColor, rect, penWidth, [on=1], [rotation=45])
%
% Draws 4 corners of a rectangle outline specified by a PTB rect of the form
% [left top right bottom]. Does this by drawing a rectangle outline and then
% drawing filled rectangles that are a bit thinner and taller / a bit
% shorter and wider in the background color, on top of the outline
% rectangle.
%
% Also rotates the placeholders by the requested rotation angle, in
% degrees.
%
% Modified from drawPlaceholders.m
%
% Rachel Denison
% March 2015


if nargin<7
    rotation = 45;
end
if nargin<6
    on = 1;
end


if on
    Screen('FramePoly', window, color, rotateCoords2(rect2Points(rect), rotation)', penWidth);
    Screen('FillPoly', window, bgColor, rotateCoords2(rect2Points(squeezeRect(rect, [0.8 1.2])), rotation)', 1); % tall rect
    Screen('FillPoly', window, bgColor, rotateCoords2(rect2Points(squeezeRect(rect, [1.2 0.8])), rotation)', 1); % wide rect
end
