function p = temporalAttentionParams

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
        p.screenRes = [1024 768];
        p.viewDist = 56;
        p.eyeTracking = 1; 
    otherwise
        error('Testing location not found in temporalAttentionParams.')
end     

p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.5;
p.goCueColor = 0.75;
p.nReps = 1;
p.nTrialsPerBlock = 64;
p.font = 'Verdana';
p.fontSize = 24;
p.showPlaceholders = 1;
p.phLineWidth = 2; % (pixels)
p.eyeRad = 1.5; % allowed fixation radius (degrees)

% Conditions
p.targetContrasts = [.64 1]; % [.64 1];
p.respInterval = [1 2]; % [1=early 2=late]
p.cueValidity = [1 -1 0]; % [1=valid -1=invalid 0=neutral]
% p.propValid = 0.67;
% p.cueValidityFactor = generatePropFactor(p.propValid);
p.cueValidityFactor = [1 1 1 2 3]; % eg. [1 1 2 3] is 50% valid, 25% invalid, 25% neutral
% p.cueValidityFactor = 1;
p.propValid = nnz(p.cueValidityFactor==1)./nnz(p.cueValidityFactor<3);
p.propNeutral = nnz(p.cueValidityFactor==3)./numel(p.cueValidityFactor);

% Timing
p.soas = [1000 1250]/1000; % [short long]
p.preCueDur = 0.75; % time between fixation onset and cue
p.cueDur = 0.2;
p.targetDur = 3/100; % 30 ms / 33 ms
p.maskSOA = 6/100; % 4/60 time between target onset and backward mask onset 
p.maskDur = 1/100; % 1/60, 3/60
p.respCueSOA = p.soas(2) + 0.5;
p.respGoSOA = 0; % 0.6 % time between resp cue onset and go onset. set to zero for no go cue.
p.iti = 0.5; % inter-trial interval (also, the duration of the feedback symbol)
p.eyeSlack = 0.12; % cushion between last fixation check and next stimulus presentation

% Images
p.imPos = [4 4];
p.imSize = [4 4]; % this is the size of the image container that holds the stim
p.targetSize = 0.5; % 0.5 sigma of gaussian / 1.5 side length of T/L / 1.5 width of triangle
p.spatialFrequency = 4; % 4
p.targetOrientation = [-5 5]; % eg. [-10 10]
p.targetPhases = 0; % eg. 0, or [0 pi/2 pi 3*pi/2]
p.TL = [0 0.5]; % [offset-for-T(=0) offset-for-L]
p.TLLineWidth = 5; % (pixels)
p.triangleBlurProp = 0.3;
p.triHWRatio = 1; % should be >= 1 so that targetSize corresponds to triangle width
p.extraOblTilt = [-5 5];
p.tiltJitter = 1;

% Staircase (implemented only for targetOrientation for now)
p.staircase = 0;
p.stairs = [.5 1 1.5 2 3 4 6 8 12];
if p.staircase
    p.targetOrientation = [0 0];
    fprintf('\nStaircase is ON\n')
end

% Task
p.task = 'targetOrientation'; % 'targetOrientation','spatialFrequency','TL'
p.targetStates = p.(p.task);
% check states
if numel(p.targetStates)~=2
    error('p.targetStates should have exactly 2 elements. Check task settings.')
end
% target rotation
switch p.task
    case 'targetOrientation'
        p.rotateTarget = 'cb'; % 'none','rotT2'= rotate T2 90 deg, 'rotT1' = rotate T1 90 deg, 'cb'= counterbalanced vert/horiz, 'cbvar'= cb with tilt variability, 'cbvark'= cb with +/-k tilt variability, 'cardobl'= cardinal + oblique orientations, 'card4'= 4 cardinal directions (triangle), 'card4wap'= like card4 but rotate the entire aperture ("card4 w/ aperture")
    case 'spatialFrequency'
        p.rotateTarget = 'random'; % random rotations
    case 'TL'
        p.rotateTarget = 'random';
        p.targetSize = 1.5;
end
%%% NB. currently, rotated targets can obscure placeholders
switch p.rotateTarget
    case {'card4','card4wap'};
        p.aperture = 'triangle'; % 'gaussian','triangle'
        p.targetSize = 2;
    otherwise
        p.aperture = 'gaussian';
end

% Masks
p.maskType = 'none'; % none, whitenoise, verticalgrating, crossedgratings, filterednoise
p.maskContrast = 1;
p.forwardMask = 0; % 1 to use forward mask, 0 for no forward mask
p.forwardMaskSOA = p.maskSOA - p.targetDur + p.maskDur; % equates ISIs between targets and masks

% Sounds
p.Fs = 44100;
% p.cueFreqs = [784 523]; % [higher G = target 1, lower C = target 2]
p.cueFreqs = [1046.5 440]; % [higher high C = target 1, lower A = target 2]
for iTone = 1:numel(p.cueFreqs)
    tone = MakeBeep(p.cueFreqs(iTone), p.cueDur, p.Fs);
    p.cueTones(iTone,:) = applyEnvelope(tone, p.Fs);
end
% 10^0.5 for every 10dB




