function y = applyEnvelope(y)
%
% function y = applyEnvelope(y)
%
% ramps a vector y on and off
%
% Rachel Denison
% March 2014

ramp = 0:.01:1; % envelope to reduce clicks
envelope = [ramp ones(1,length(y)-length(ramp)*2) 1-ramp];
y = y.*envelope;
