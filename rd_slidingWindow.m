function [ySmoothed, steps] = rd_slidingWindow(x, y, winSize, xrange)
%
% function [ySmoothed, steps] = rd_slidingWindow(x, y, winSize, xrange)
%
% example inputs:
% x = targetOrientDiff{1,1};
% y = errors{1,1};
% winSize = 15; % should be odd
% xrange = [-90 90];

steps = xrange(1)+floor(winSize/2):xrange(2)-floor(winSize/2);
ySmoothed = nan(1,numel(steps));
for i = 1:numel(steps)
    win = steps(i)-floor(winSize/2):steps(i)+floor(winSize/2);
    w = [];
    for j = 1:winSize
        w(:,j) = x==win(j);
    end
    wAll = any(w,2);
    ySmoothed(i) = mean(y(wAll));
end

% figure
% hold on
% plot(x, y, '.')
% plot(steps, ySmoothed, 'r', 'LineWidth', 2)