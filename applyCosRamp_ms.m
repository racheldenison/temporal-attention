function [s] = applyCosRamp_ms(sig,durRamp,sf)
%APPLYCOSRAMP_MS This function applies a hanning window cosine ramp at 
% the first and last x ms of the signal.
%
%   INPUTS:
%   sig:- signal to be adjusted.
%   durRamp:- duration of signal to be ramped at onset/offset in ms.
%   sf:- sampling frequency.
%
%   OUTPUTS:
%   s:- ramped signal.
%
% Created by SML Mar 2015

assert(isvector(sig),'Please input a vector for sig.')

adj = round(durRamp/1000 * sf); % #samples for each ramp
win = hanning(adj*2); % Hanning window cosine ramp

% If necessary, rotate win vector to match signal:
[temp,c_sig]=size(sig);
[temp,c_win]=size(win);
if c_sig ~= c_win
   win = win'; 
end

% Onset and offset components:
win_on = win(1:adj);
win_off = win(adj+1:end);

% Apply ramp to signal:
s = sig;
s(1:adj) = s(1:adj) .* win_on;
s(end-adj+1:end) = s(end-adj+1:end) .* win_off;

end