function p = temporalAttentionParams3Targets

p.testingLocation = 'CarrascoL1'; % 'CarrascoL1','laptop','desk'

switch p.testingLocation
    case {'laptop','desk'}
        p.keyNames = {'1!','2@'};
        p.refRate = 1/60;
        p.screenSize = [9 13]; % (in)
        p.screenRes = [900 1440];
        p.viewDist = 21; % (in)
        p.eyeTracking = 0;
    case 'CarrascoL1'
        p.keyNames = {'1!','2@'};
        p.refRate = 1/100;
        p.screenSize = [40 30];
        p.screenRes = [1280 960]; % [1024 768]
        p.viewDist = 57; % 56
        p.eyeTracking = 1; 
    otherwise
        error('Testing location not found in temporalAttentionParams.')
end

p.sessionStartBlocks = [6 11];
p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.5;
p.goCueColor = 0.75;
p.nReps = 1;
p.nTrialsPerBlock = 64; % 64
p.font = 'Verdana';
p.fontSize = 24;
p.showPlaceholders = 1;
p.phLineWidth = 2; % (pixels)
p.eyeRad = 1.5; % allowed fixation radius (degrees)    

% Condition
p.targetContrasts = 1; % [.64 1] [.16 .64];
p.respInterval = [1 2 3]; % [1=T1 2=T2 3=T3]
p.cueValidity = [1 -1 0]; % [1=valid -1=invalid 0=neutral]
% p.propValid = 0.67;
% p.cueValidityFactor = generatePropFactor(p.propValid);
p.cueValidityFactor = [1 1 1 2 3]; % eg. [1 1 2 3] is 50% valid, 25% invalid, 25% neutral
% p.cueValidityFactor = 1;
p.propValid = nnz(p.cueValidityFactor==1)./nnz(p.cueValidityFactor<3);
p.propNeutral = nnz(p.cueValidityFactor==3)./numel(p.cueValidityFactor);

% Timing
p.soas = [1000 1250 1500]/1000; % [T1 T2 T3]
p.preCueDur = 0.75; % time between fixation onset and cue
p.cueDur = 0.2;
p.targetDur = 3/100; % 30 ms / 33 ms
p.maskSOA = p.soas(2) - p.soas(1);%4/100; % 4/60 time between target onset and backward mask onset %%%% 1/100 to match other runs vp_cbD6
p.maskDur = p.targetDur; % 1/60, 3/60
% p.forwardMaskSOA = p.maskSOA - p.targetDur + p.maskDur; % equates ISIs between targets and masks
p.forwardMaskSOA = p.maskSOA; 
p.respCueSOA = p.soas(end) + 0.5;
p.respGoSOA = 1.5; % 0.6 % 1.0 % time between resp cue onset and go onset. set to zero for no go cue.
p.iti = 0.5; % inter-trial interval (also, the duration of the feedback symbol)
p.eyeSlack = 0.12; % cushion between last fixation check and next stimulus presentation

% Images
p.imPos = [4 4];
p.imSize = [4 4]; % this is the size of the image container that holds the stim
p.targetSize = 0.5; % 0.5 sigma of gaussian / 1.5 side length of T/L / 1.5 width of triangle
p.spatialFrequency = 4; % 4
p.targetOrientation = [-4 4]; % eg. [-10 10]
p.targetPhases = 0; % eg. 0, or [0 pi/2 pi 3*pi/2]

% Staircase (implemented only for targetOrientation for now)
p.staircase = 0;
p.stairs = [.5 1 1.5 2 3 4 6 8 12]; % 12
if p.staircase
    p.targetOrientation = [0 0];
    fprintf('\nStaircase is ON\n')
end

% Task
p.task = 'targetOrientation'; % 'targetOrientation','spatialFrequency'
p.targetStates = p.(p.task);
% check states
if numel(p.targetStates)~=2
    error('p.targetStates should have exactly 2 elements. Check task settings.')
end
% target rotation
switch p.task
    case 'targetOrientation'
        p.rotateTarget = 'cb'; % 'none','random','cb'= counterbalanced vert/horiz,'vh'= vert/horiz not fully counterbalanced
    case 'spatialFrequency'
        p.rotateTarget = 'random'; % random rotations
end
%%% NB. currently, rotated targets can obscure placeholders
p.aperture = 'gaussian';

% Masks
p.maskType = 'none'; % none, whitenoise, verticalgrating, crossedgratings, filterednoise, bullseye, pseudotarget, h&vgratings
p.maskContrast = 1;
p.maskSFBand = [1/1.3 1.3]*p.spatialFrequency;
p.forwardMask = [0 0 0]; % T1, T2, T3     1 to use forward mask, 0 for no forward mask
p.backwardMask = [0 0 0]; % T1, T2, T3    1 to use backward mask, 0 for no backward mask

% Sounds
p.Fs = 44100;
p.cueFreqs = [784 523 330]; % [higher G5 = target 1, lower C5 = target 2, lowest E4 = target 3]
% p.cueFreqs = [1046.5 440]; % [higher high C = target 1, lower A = target 2]
for iTone = 1:numel(p.cueFreqs)
    tone = MakeBeep(p.cueFreqs(iTone), p.cueDur, p.Fs);
    p.cueTones(iTone,:) = applyEnvelope(tone, p.Fs);
end
p.playTargetSound = 0;
p.targetSoundDur = 0.03;
p.targetSoundFreq = [200 10000];
p.targetSound = makeNoiseBurst(0, p.targetSoundDur, 5,...
    p.targetSoundFreq(1), p.targetSoundFreq(2), p.Fs, 0)';
% p.targetSound = applyEnvelope(MakeBeep(p.targetSoundFreq, p.targetSoundDur, p.Fs), p.Fs);
% 10^0.5 for every 10dB




