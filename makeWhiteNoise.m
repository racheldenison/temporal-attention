function [s_all,noise_all] = makeWhiteNoise(nSounds,selFilter,lowF,highF,dur,sf,plotyn)
%MAKEWHITENOISE This function generates a single channel white noise
%stimulus with filtering of your choice (see below).
%
%   INPUTS:
%   nSounds:- number of white noise stimuli
%   filter:- (0) None, (1) Lowpass, (2) Highpass, (3) Bandpass, (4) Notch.
%   lowF:- lowest frequency
%   highF:- highest frequency
%   dur:- duration in seconds of sound stimulus
%   sf:- sampling frequency
%
%   OUTPUTS:
%   noise_all:- the pre-filtered noise
%   s_all:- the filtered sound (1 channel)
%
%   HISTORY:
%   1. Code sourced from: http://www.h6.dion.ne.jp/~fff/old/technique/auditory/matlab.html
%   2. Adapted by SML Oct 2014


%---------------
% Set Parameters
%---------------

% Check if filter option selected is valid:
assert(ismember(selFilter,0:4),'Select a filter (0-4).')

% Set defaults:
if nargin < 7
    plotyn = 0;
    if nargin < 6
        sf = 48000;
        if nargin < 5
            dur = 1;
            if nargin < 4
                highF = 16000;
                if nargin < 3
                    lowF = 300;
                    if nargin < 2
                        selFilter = 0;
                        if nargin < 1
                            nSounds = 1;
                        end
                    end
                end
            end
        end
    end
end

% General variables:
nSamples = sf * dur;  % number of samples
nf = sf / 2; % nyquist frequency
nh = nSamples / 2;  % half number of samples

% Filter variables:
lp = lowF * dur; % lowF point in frequency domain
hp = highF * dur; % highF point in frequency domain

%--------------
% Design Filter
%--------------

switch selFilter
    case 0
        fType = 'NONE';
    case 1
        fType = 'LOWPASS';
        filter = zeros(1, nSamples);           % initializaiton by 0
        filter(1, 1 : lp) = 1;          % filter design in real number
        filter(1, nSamples - lp : nSamples) = 1;      % filter design in imaginary number
    case 2
        fType = 'HIGHPASS';
        filter = ones(1, nSamples);            % initializaiton by 1
        filter(1, 1 : hp) = 0;          % filter design in real number
        filter(1, nSamples - hp : nSamples) = 0;      % filter design in imaginary number
    case 3
        fType = 'BANDPASS';
        filter = zeros(1, nSamples);           % initializaiton by 0
        filter(1, lp : hp) = 1;         % filter design in real number
        filter(1, nSamples - hp : nSamples - lp) = 1; % filter design in imaginary number
    case 4
        fType = 'NOTCH';
        filter = ones(1, nSamples);
        filter(1, lp : hp) = 0;
        filter(1, nSamples - hp : nSamples - lp) = 0;
end

% Display the filter selected in command window:
% fprintf('\n\n Filter used: '); disp(fType)

%----------------------
% Make Noise and Filter
%----------------------

noise_all = zeros(nSamples,nSounds); % Preallocate
s_all = zeros(nSamples,nSounds); % Preallocate
rand('state',sum(100 * clock));  % initialize random seed

for ii = 1:nSounds
    
    % Noise:
    noise = randn(1, nSamples);             % Gausian noise
    noise = noise / max(abs(noise)); % -1 to 1 normalization
    noise_all(:,ii) = noise';
    
    % Filter:
    if selFilter ~= 0
        s = fft(noise);                  % FFT
        s = s .* filter;                 % filtering
        s = ifft(s);                     % inverse FFT
        s = real(s);
    else
        s = noise;
    end
    s_all(:,ii) = s';
    
end


%----------------------------
% Plot Example Filtered Sound
%----------------------------
% (Will plot the last generated stimulus)

if plotyn == true
    
    % Plot sounds:
    x = linspace(0, dur, nSamples);
    subplot(2,2,1); plot(x, noise); xlabel('time (s)'); title('sound: noise');
    subplot(2,2,2); plot(x, s); xlabel('time (s)'); title(['sound: filtered noise, ' fType]);
    
    % Plot Fourier spectrums:
    x = linspace(0, nf, nh);
    t = fft(noise);
    t = t .* conj(t);
    subplot(2,2,3); semilogy(x, t(1,1:nh) ./ max(t));  xlabel('frequency (Hz)'); title('spectrum: noise');
    t = fft(s);
    t = t .* conj(t);
    subplot(2,2,4); semilogy(x, t(1,1:nh) ./ max(t));  xlabel('frequency (Hz)');  title(['spectrum: filtered noise, ' fType]);
    
    figure(1);
    
end

end