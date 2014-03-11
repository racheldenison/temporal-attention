function p = temporalAttentionParams

p.testingLocation = 'laptop'; % 'CarrascoL1','laptop'

switch p.testingLocation
    case 'laptop'
        p.keyNames = {'1!','2@'};
        p.refRate = 1/60;
        p.screenSize = [9 13]; % (in)
        p.screenRes = [900 1440];
        p.viewDist = 21; % (in)
    case 'CarrascoL1'
        p.keyNames = {'1!','2@'};
        p.refRate = 1/60;
        p.screenSize = [40 30];
        p.screenRes = [1024 768];
        p.viewDist = 57;
    otherwise
        error('Testing location not found in temporalAttentionParams.')
end     

p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.5;
p.nReps = 1;
p.font = 'Verdana';

% Conditions
p.targetContrasts = [0.16 0.32];
p.respInterval = [1 2]; % [1=early 2=late]
p.cueValidity = [1 -1]; % [1=valid -1=invalid]
p.propValid = 0.75;
p.cueValidityFactor = generatePropFactor(p.propValid);

% Timing
p.soas = [1000 1433]/1000; % [short long]
p.cueDur = 0.2;
p.targetDur = 2/60; % 33 ms
p.respCueSOA = p.soas(2) + 0.5;
p.iti = 1; % inter-trial interval

% Images
p.fixSize = 0.1;
p.imPos = [4 4];
p.imSize = [6 6]; % this is the size of the image container that holds the stim
p.targetSize = 0.5; % sigma of gaussian
p.spatialFrequency = 4;
p.targetOrientations = [-15 15];

% Sounds
p.Fs = 44100;
p.cueFreqs = [784 523]; % [higher G = target 1, lower C = target 2]
t = 0:1/p.Fs:1-1/p.Fs; % make a 1-s vector with sampling freq Fs
t = t(1:p.Fs*p.cueDur); % truncate t at cueDur
for iTone = 1:numel(p.cueFreqs)
    p.cueTones(iTone,:) = sin(2*pi*p.cueFreqs(iTone)*t);
end
% soundsc(p.cueTones(1,:), p.Fs)
% 10^0.5 for every 10dB
%
% Old working example:
% v = 1:1000;
% envelope = [0:.05:1 ones(1,1000-42) 1:-.05:0];
% startSound = (0.25 * sin(3*pi*v/30)).*envelope;





