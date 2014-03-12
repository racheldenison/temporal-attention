function y = applyEnvelope(y, fs, rampDur_ms)
%
% function y = applyEnvelope(y, fs, [rampDur_ms])
%
% ramps a vector y on and off
% fs is the audio sampling frequency
% rampDur_ms is the duration of each ramp in ms (optional: defaults to 10)
%
% Rachel Denison
% March 2014
%
% Cosine ramp section from Elise Piazza

%% Linear ramp
% ramp = 0:.01:1; % envelope to reduce clicks
% envelope = [ramp ones(1,length(y)-length(ramp)*2) 1-ramp];
% y = y.*envelope;

%% Cosine ramp
if nargin<3
    rampDur_ms = 10; % ramp duration (beginning and end) cosinus: 10 ms
end

% computation of cosine ramps %
rampDur_samples = round(fs*(rampDur_ms/1000)); % duration  of  ramps in samples
Freq_ramp = 1/(2*(rampDur_ms/1000)); % frequence of the cosine gate
temps = 1:rampDur_samples; % generate ramp vector
onset  = (1 + sin(2*pi*Freq_ramp*temps./fs + (-pi/2)))/2;
offset = (1 + sin(2*pi*Freq_ramp*temps./fs + (pi/2)))/2;

% define env sound - here I just add the onset and offset ramps to the
% sound by multiplying the onset or offset vector with my sound vector (for
% more complicated amplitude modulations: just define another onset or
% offset vector)
y(1:length(onset)) =  y(1:length(onset)) .* onset;
y(end-length(offset)+1:end) =  y(end-length(offset)+1:end) .* offset;
        