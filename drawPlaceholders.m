function drawPlaceholders(window, color, bgColor, rect, penWidth, on)

% function drawPlaceholders(window, color, bgColor, rect, penWidth, [on=1])
%
% Draws 4 corners of a square outline specified by a PTB rect of the form
% [left top right bottom]
%
% Rachel Denison
% March 2014

if nargin<6
    on = 1;
end

if on
    Screen('FrameRect', window, color, rect, penWidth);
    Screen('FillRect', window, bgColor, squeezeRect(rect, [0.8 1.2])) % tall rect
    Screen('FillRect', window, bgColor, squeezeRect(rect, [1.2 0.8])) % wide rect
end