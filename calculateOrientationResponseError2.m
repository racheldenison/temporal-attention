function [err errAbs] = calculateOrientationResponseError2(rot, resp)
%
% identical to calculateOrientationResponseError, but works with vectors

err = resp - rot;
err = mod((err + 180), 360) - 180;
err(err < -90) = err(err < -90) + 180;
err(err > 90) = err(err > 90) - 180;
errAbs = abs(err);
