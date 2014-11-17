function [err errAbs] = calculateOrientationResponseError(rot, resp)

% err = resp - rot;
% 
% % fixes case of obtuse angle on the 0-180 range
% if abs(err) > 90
%     err = resp - rot - 180;
% end
% 
% % fixes case of crossing 0/180
% if abs(err) > 90
%     err = 180 - rot + resp;
% end
% 
% errAbs = abs(err);

% simpler
err = resp - rot;
err = mod((err + 180), 360) - 180;
if err < -90
    err = err + 180;
elseif err > 90
    err = err - 180;
end
errAbs = abs(err);
