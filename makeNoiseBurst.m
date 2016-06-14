function [click] = makeNoiseBurst(playSoundYN,dur,cosWin,lowF,highF,sf,durPadding)
% MAKENOISEBURST will create a band-passed click stimulus.
%
% INPUTS:
% playSoundYN = Do you want this script to play the generated sound?
% dur = duation of click stimulus (in seconds)
% cosWin = rise-fall time (in ms)
% lowF = low frequency cut off
% highF = ligh frequency cut off
% sf = sampling frequency of stimulus
% durPadding = The duration of leading/trailing silence (it's important to
% pad the auditory buffer)
%
% OUTPUTS:
% s = click sound
%
% Written by SML June 2016

if nargin < 7
    durPadding = 0.005;
    if nargin < 6
        sf = 48000;
        if nargin < 5
            highF = 10000;
            if nargin < 4
                lowF = 200;
                if nargin < 3
                    cosWin = 5;
                    if nargin < 2
                        dur = 0.02;
                        if nargin < 1
                            playSoundYN = 1;
                        end
                    end
                end
            end
        end
    end
end

% Create a bandpass noise stimulus:
click = makeWhiteNoise(1,3,lowF,highF,dur,sf,0);

% Apply onset-offset ramp for smooth sound:
click = applyCosRamp_ms(click,cosWin,sf);

% Apply the zero padding:
click = [zeros(round(durPadding*sf),1); click; zeros(round(durPadding*sf),1)];

% Play the sound if option selected:
% (note that MATLAB will continue executing script once the sound has been
% begins playing)
if playSoundYN == 1
    player = audioplayer(click,sf);
    play(player);
    WaitSecs(0.5);
end

end