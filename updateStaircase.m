function [stairIdx lastFewAcc] = updateStaircase(stairs, stairIdx, lastFewAcc, correct)

% [stairIdx lastFewAcc] = updateStaircase(stairs, stairIdx, lastFewAcc, correct)
%
% Implements a 3 up, 1 down staircase
% Output args are required and must be named the same as their
% corresponding input args.
% Assumes the stairs are ordered hard to easy (as when finding a threshold)
%
% Rachel Denison
% 2014 April 15

% keep track of last 3 trials
lastFewAcc = [lastFewAcc correct];
if numel(lastFewAcc)>3
    lastFewAcc = lastFewAcc(end-2:end);
end
if sum(lastFewAcc)==3
    stairIdx = stairIdx-1; % make it harder
elseif correct==0
    stairIdx = stairIdx+1; % make it easier
end
if stairIdx > numel(stairs)
    stairIdx = numel(stairs);
elseif stairIdx < 1
    stairIdx = 1;
end
