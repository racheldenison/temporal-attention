function results = rd_temporalAttention(subjectID, varargin)

%% Setup
if nargin<1
    subjectID = 'test';
end

saveData = 1;
saveFigs = 1;
plotTimingFigs = 1;
saveTimingFigs = 1;

if nargin<2
    workspaceFile = [];
end
if nargin>1
    workspaceOption = varargin{1};
    if isempty(workspaceOption)
        workspaceFile = [];
    elseif ischar(workspaceOption) && strcmp(workspaceOption,'w')
        % workspaceFile = [pathToExpt('data') '/pilot/ad/adPilot_cb_tc64-100_soa1000-1250_run01_WORKSPACE.mat'];
        workspaceFile = sprintf('data/%s_WORKSPACE.mat', subjectID);
    else
        error('to load from workspace file, input "w" as second argument')
    end
end
if nargin>2
    p0 = varargin{2};
end
if nargin>3
    error('wrong number of input arguments')
end

p = temporalAttentionParams;

% note, p fields may have dependencies: check that all fields are set!
if exist('p0','var')
    fNames = fields(p0);
    for iF = 1:numel(fNames)
        fName = fNames{iF};
        p.(fName) = p0.(fName);
    end
end
    
if strcmp(subjectID,'test')
    p.eyeTracking = 0;
end

slack = p.refRate/2;

rad = round(ang2pix(p.eyeRad, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

%% Display key settings to the experimenter
fprintf('\nExperiment settings:\n')
fprintf('tilt = [%1.1f %1.1f]\n', p.targetOrientation(1), p.targetOrientation(2))
fprintf('soa = %d ms\n', round(1000*(p.soas(2)-p.soas(1))))
fprintf('respGoSOA = %d ms\n\n', 1000*p.respGoSOA)

ok = input('Settings ok? [n if not]','s');
if strcmp(ok,'n')
    error('Ok, check parameters')
end

%% Eye data i/o
eyeDataDir = 'eyedata';
eyeFile = sprintf('%s%s', subjectID([1:2 end-1:end]), datestr(now, 'mmdd'));

% Check to see if this eye file already exists
existingEyeFile = dir(sprintf('%s/%s.edf', eyeDataDir, eyeFile));
if ~isempty(existingEyeFile) && p.eyeTracking
%     error('eye file already exists! please choose another name.')
end

%% Check for existing data file
% note the real data file is made in the analyze function
dataFile = sprintf('%s_TemporalAttention_T1T2all_%s.mat', subjectID, datestr(now, 'yyyymmdd'));
existingDataFile = dir(sprintf('data/%s', dataFile));
if ~isempty(existingDataFile) && ~strcmp(subjectID, 'test') && ~strcmp(subjectID, 'testy')
    error('data file already exists!')
end

%% Keyboard
% Find keyboard device number
devNum = findKeyboardDevNumsAtLocationNYU(p.testingLocation);
if isempty(devNum)
    %error('Could not find Keypad!')
%     devNum = -1;
end
devNum = -1;

%% Sound
% Perform basic initialization of the sound driver
InitializePsychSound(1); % 1 for precise timing

% Open audio device for low-latency output
reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel

%% Screen
% Set up screen
screenNumber = max(Screen('Screens'));

% Check screen resolution and refresh rate - if it's not set correctly to 
% begin with, the color might be off
scr = Screen('Resolution', screenNumber);
if ~all([scr.width scr.height scr.hz] == [p.screenRes round(1/p.refRate)]) && ~strcmp(subjectID,'test')
    error('Screen resolution and/or refresh rate has not been set correctly by the experimenter!')
end

% % Set resolution and refresh rate
% if strcmp(p.testingLocation, 'CarrascoL1')
%     Screen('Resolution',screenNumber, p.screenRes(1), p.screenRes(2), round(1/p.refRate));
% end

% Set up window
[window rect] = Screen('OpenWindow', screenNumber);
% [window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]);
white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
[cx cy] = RectCenter(rect);
Screen('TextSize', window, p.fontSize);
Screen('TextColor', window, white);
Screen('TextFont', window, p.font);

% Check screen size
[sw, sh] = Screen('WindowSize', window); % height and width of screen (px)
if ~all([sw sh] == p.screenRes) && ~strcmp(subjectID,'test')
    error('Screen resolution is different from requested!')
end

% Check refresh rate
flipInterval = Screen('GetFlipInterval', window); % frame duration (s)
if abs(flipInterval - p.refRate) > 0.001 && ~strcmp(subjectID,'test')
    error('Refresh rate is different from requested!')
end

% Check font
if ~strcmp(p.font, Screen('TextFont', window))
    error('Font was not set to requested: %s', p.font)
end

% Load calibration file
switch p.testingLocation
    case 'CarrascoL1'
%         load('../../Displays/0001_james_TrinitonG520_1280x960_57cm_Input1_140129.mat');
%         table = repmat(calib.table,1,3);
        load('../../Displays/Carrasco_L1_SonyGDM5402_sRGB_calibration_02292016.mat');
        table = CLUT;
        Screen('LoadNormalizedGammaTable', window, table);
        % check gamma table
        gammatable = Screen('ReadNormalizedGammaTable', window);
        if nnz(abs(gammatable-table)>0.0001)
            error('Gamma table not loaded correctly! Perhaps set screen res and retry.')
        end
    otherwise
        fprintf('\nNot loading gamma table ...\n')
end

%% Make stimuli
% Calculate stimulus dimensions (px) and position
imPos = round(ang2pix(p.imPos, p.screenSize(1), p.screenRes(1), p.viewDist, 'radial')); % from screen center
targetSize = round(ang2pix(p.targetSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
pixelsPerDegree = round(ang2pix(1, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));

% Make target gratings
for iP = 1:numel(p.targetPhases)
    targetPhase = p.targetPhases(iP);
    for iC = 1:numel(p.targetContrasts)
        targetContrast = p.targetContrasts(iC);
        for iT = 1:numel(p.targetStates)
            targetState = p.targetStates(iT);
            
            switch p.task
                case 'TL'
                    targetColor = targetContrast*(1-p.backgroundColor) + p.backgroundColor;
                    t = makeTLImage(targetSize, p.TLLineWidth, targetState, targetColor, p.backgroundColor);
                otherwise
                    switch p.task
                        case 'targetOrientation'
                            targetOrientation = targetState;
                            spatialFrequency = p.spatialFrequency;
                        case 'spatialFrequency'
                            spatialFrequency = targetState;
                            targetOrientation = p.targetOrientation;
                        otherwise
                            error('p.task not recognized')
                    end
                    
                    % overwrite targetOrientation if we are also rotating the
                    % aperture. will take care of the target tilt via
                    % targetRotations.
                    if strcmp(p.rotateTarget, 'card4wap')
                        targetOrientation = 0;
                    end
                    
                    % make big grating
                    t = buildColorGrating(pixelsPerDegree, p.imSize, ...
                        spatialFrequency, targetOrientation, targetPhase, targetContrast, 0, 'bw');
            end
            
            % mask with an aperture (eg. 2d gaussian)
            switch p.aperture
                case 'gaussian'
                    target{iC,iT,iP} = maskWithGaussian(t, size(t,1), targetSize);
                case 'triangle'
                    blurSize = round(targetSize*p.triangleBlurProp);
                    target{iC,iT,iP} = maskWithTriangle(t, ...
                        round(size(t,2)/2), round(size(t,1)/2), ...
                        targetSize, targetSize*p.triHWRatio, blurSize);
                otherwise
                    error('p.aperture not recognized')
            end
        end
    end
end

% Make mask
switch p.maskType
    case 'none'
%         mask = ones(size(t))*p.backgroundColor;
        mask = ones(size(t)); % white so that it will be obvious if it is being presented when it shouldn't be
    case 'whitenoise'
        mask = (rand(size(t))-0.5)*p.maskContrast + 0.5;
    case 'verticalgrating'
        m = buildColorGrating(pixelsPerDegree, p.imSize, ...
            p.spatialFrequency, 0, 0, p.maskContrast, 0, 'bw');
        mask = maskWithGaussian(m, size(m,1), targetSize);
    case 'crossedgratings'
        for i=1:numel(p.targetOrientation)
            m(:,:,i) = buildColorGrating(pixelsPerDegree, p.imSize, ...
                p.spatialFrequency, p.targetOrientation(i), 0, p.maskContrast, 0, 'bw');
        end
        m = sum(m,3)./numel(p.targetOrientation);
        mask = maskWithGaussian(m, size(m,1), targetSize);
    case 'hvgratings'
        hv = [0 90];
        for i=1:numel(hv)
            m(:,:,i) = buildColorGrating(pixelsPerDegree, p.imSize, ...
                p.spatialFrequency, hv(i), 0, p.maskContrast, 0, 'bw');
        end
        m = sum(m,3)./numel(hv);
        mask = maskWithGaussian(m, size(m,1), targetSize);
    case 'filterednoise'
        idx = 1;
        while idx <= 100
%             masktemp = makeFilteredNoise(p.imSize(1)/1.3, p.maskContrast, ...
%                 0, 180, p.spatialFrequency, 2, pixelsPerDegree, 1);
            masktemp = makeFilteredNoise2(p.imSize(1), p.maskContrast, ...
                0, 180, p.maskSFBand(1), p.maskSFBand(2), pixelsPerDegree, 0);
            masktemp = maskWithGaussian(masktemp, size(masktemp,1), targetSize);
            if mean(masktemp(:))<0.51 && mean(masktemp(:))>0.49
                % match rms contrast
                for i = 1:5 % a few loops of scaling and truncation
                    maskRMSFactor = std(target{1,1,1}(:))/std(masktemp(:))
                    masktemp = (masktemp-.5)*maskRMSFactor + .5;
                    masktemp(masktemp>1) = 1;
                    masktemp(masktemp<0) = 0;
                end
                mask{idx} = masktemp;
                idx = idx+1;
            end
        end
    case 'bullseye'
       m = buildBullseye(p.spatialFrequency,p.imSize,pixelsPerDegree,p.maskContrast);
       mask = maskWithGaussian(m, size(m,1), targetSize);
    otherwise
        error('maskType not recognized')
end

%% Make textures
for iP = 1:numel(p.targetPhases)
    for iC = 1:numel(p.targetContrasts)
        for iT = 1:numel(p.targetStates)
            targetTex(iC,iT,iP) = Screen('MakeTexture', window, target{iC,iT,iP}*white);
        end
    end
end

if iscell(mask)
    for i=1:numel(mask)
        maskTexs(i) = Screen('MakeTexture', window, mask{i}*white);
    end
else
    maskTexs = Screen('MakeTexture', window, mask*white);
end

% Make the rects for placing the images
imSize = size(target{1,1},1);
imRect = CenterRectOnPoint([0 0 imSize imSize], cx+imPos(1), cy+imPos(2));
phRect = imRect + [-1 -1 1 1]*p.phLineWidth;
phRect = phRect + [-1 -1 1 1]*round(ang2pix(0.25, p.screenSize(1), p.screenRes(1), p.viewDist, 'central')); % expand placeholders by this many pixels so that they are not obscured by image rotations

%% Generate trials
% Construct trials matrix
% a "state" is one of two states the target can take, eg left vs. right
% orientation, or low vs. high SF
trials_headers = {'targetContrast','respInterval','cueValidity',...
    'target1State','target2State','cuedInterval','respTargetState',...
    'rt','responseKey','response','correct'};

% make sure column indices match trials headers
targetContrastIdx = strcmp(trials_headers,'targetContrast');
respIntervalIdx = strcmp(trials_headers,'respInterval');
cueValidityIdx = strcmp(trials_headers,'cueValidity');
target1StateIdx = strcmp(trials_headers,'target1State');
target2StateIdx = strcmp(trials_headers,'target2State');
cuedIntervalIdx = strcmp(trials_headers,'cuedInterval');
respTargetOrientIdx = strcmp(trials_headers,'respTargetState');
rtIdx = strcmp(trials_headers,'rt');
responseKeyIdx = strcmp(trials_headers,'responseKey');
responseIdx = strcmp(trials_headers,'response');
correctIdx = strcmp(trials_headers,'correct');

% full factorial design
trials = fullfact([numel(p.targetContrasts) ...
    numel(p.respInterval) ...
    numel(p.cueValidityFactor) ...
    numel(p.targetStates) ...
    numel(p.targetStates)]);
nBaseCond = size(trials,2);

% generate target rotations
% these are applied over and above any initial orientation
nTrials0 = size(trials,1);
switch p.rotateTarget
    case 'none'
        targetRotations = zeros(nTrials0,2);
    case 'random'
        targetRotations = 360*rand(nTrials0,2);
    case 'rotT1'
        targetRotations = [90*ones(nTrials0,1) zeros(nTrials0,1)];
    case 'rotT2'
        targetRotations = [zeros(nTrials0,1) 90*ones(nTrials0,1)];
    case 'cb'
        % if we want to counterbalance vertical and horizontal rotations 
        % across trials, we need to have 4x the original number of trials
        targetRotations(:,1) = repmat([zeros(nTrials0,1); 90*ones(nTrials0,1)],2,1);
        targetRotations(:,2) = [zeros(nTrials0*2,1); 90*ones(nTrials0*2,1)];
        trials = repmat(trials,4,1);
    case 'cbvar'
        % same as cb, with the added twist that we will add a little bit of
        % variability about the tilt
        targetRotations(:,1) = repmat([zeros(nTrials0,1); 90*ones(nTrials0,1)],2,1);
        targetRotations(:,2) = [zeros(nTrials0*2,1); 90*ones(nTrials0*2,1)];
        trials = repmat(trials,4,1);
        % now add the extra variability around the base tilt
        baseTilt = min(abs(p.targetStates));
        tiltJitterRange = [-0.75*baseTilt .75*baseTilt];
        tiltJitters = tiltJitterRange(1) + diff(tiltJitterRange).*rand(size(targetRotations));
        targetRotations = targetRotations + tiltJitters;
    case 'cbvark'
        % same as cbvar, but we'll only use tilt, tilt+k, and tilt-k
        targetRotations(:,1) = repmat([zeros(nTrials0,1); 90*ones(nTrials0,1)],2,1);
        targetRotations(:,2) = [zeros(nTrials0*2,1); 90*ones(nTrials0*2,1)];
        trials = repmat(trials,4,1);
        % use the strategy from cardobl to get roughly equal numbers of
        % each tilt jitter pair. require that T1 and T2 never have
        % identical tilts
        tiltJitters0 = generatePseudorandomPairs(3, size(targetRotations,1), 1);
        tiltJitters = zeros(size(tiltJitters0));
        tiltJitters(tiltJitters0==1) = -p.tiltJitter;
        tiltJitters(tiltJitters0==3) = p.tiltJitter;
        targetRotations = targetRotations + tiltJitters;
    case 'cardobl'
        % if we wanted to counterbalance across all possible rotation pairs
        % in the early and late trials, we would need to have 12x the
        % original number of trials (length of rotpairs, below). This seems
        % impractical. so let's just try to have roughly equal numbers of 
        % each rot pair
        targetRotations = generatePseudorandomPairs(4, nTrials0, 1);
        targetRotations = targetRotations.*45 - 45;
        % add extra tilt to oblique orientations
        isobl = mod(targetRotations,90)~=0; % find the obliques
        tilts = trials(:,target1StateIdx | target2StateIdx);
        targetRotations = targetRotations + p.extraOblTilt(tilts).*isobl;
    case 'card4'
        % same strategy as in cardobl: dp not counterbalance, but try to
        % have the same number of each cardinal direction pair
        targetRotations = generatePseudorandomPairs(4, nTrials0, 1);
        targetRotations = targetRotations.*90 - 90;
    case 'card4wap'
        % same strategy as card4. however, all gratings will be zero tilt
        % originally, so add extra rotations as in cardobl.
        targetRotations = generatePseudorandomPairs(4, nTrials0, 1);
        targetRotations = targetRotations.*90 - 90;
        % add rotations to all targets
        tilts = trials(:,target1StateIdx | target2StateIdx);
        targetRotations = targetRotations + p.targetOrientation(tilts);
    otherwise
        error('p.rotateTarget not recognized')
end

% generate target phases (indices for the target textures)
targetPhases = generatePseudorandomPairs(numel(p.targetPhases), size(trials,1), 1);

% set cue validity condition according to the cueValidityFactor (potentially unequal proportion of trials in each condition)
cueValidityTrials = trials(:,cueValidityIdx);
for i = 1:numel(p.cueValidityFactor)   
    trials(cueValidityTrials==i,cueValidityIdx) = p.cueValidityFactor(i);
end

% repeat trials matrix according to nReps of all conditions
trials = repmat(trials, p.nReps, 1);
targetRotations = repmat(targetRotations, p.nReps, 1);
targetPhases = repmat(targetPhases, p.nReps, 1);
nTrials = size(trials,1);

% show trials and blocks
fprintf('\n%s\n\n%d trials, %1.2f blocks\n\n', datestr(now), nTrials, nTrials/p.nTrialsPerBlock)

% Choose order of trial presentation
% trialOrder = randperm(nTrials);
% This randomizes trial order within reps, but not across reps. So we will
% present every trial in one rep before moving on to the next rep. In this
% way the experiment is blocked into complete reps (though the number of
% trials per block may be different then the number of trials in a rep).
nt = nTrials/p.nReps;
trialSet = ones(nt,1)*(1:p.nReps);
repOrders = [];
for i=1:p.nReps
    repOrders = [repOrders randperm(nt)'];
end
trialOrder0 = repOrders + (trialSet-1).*nt;
trialOrder = trialOrder0(:);

%% Eyetracker
if p.eyeTracking    
    % Initialize eye tracker
    [el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
    if exitFlag
        return
    end
    
    % Write subject ID into the edf file
    Eyelink('message', 'BEGIN DESCRIPTIONS');
    Eyelink('message', 'Subject code: %s', subjectID);
    Eyelink('message', 'END DESCRIPTIONS');
    
    % No sounds indicating success of calibration
%     el.targetbeep = false;
%     el.calibration_failed_beep = [0 0 0];
%     el.calibration_success_beep = [0 0 0];
    el.drift_correction_target_beep = [0 0 0];
    el.drift_correction_failed_beep = [0 0 0];
    el.drift_correction_success_beep = [0 0 0];
    
    % Accept input from all keyboards
    el.devicenumber = -1; %see KbCheck for details of this value
    
    % Update with custom settings
    EyelinkUpdateDefaults(el);

    % Calibrate eye tracker
    [cal exitFlag] = rd_eyeLink('calibrate', window, el);
    if exitFlag
        return
    end
end

%% Present trials
% Show fixation and wait for a button press
Screen('FillRect', window, white*p.backgroundColor);
DrawFormattedText(window, 'x', 'center', 'center', white);
drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
Screen('Flip', window);
KbWait(devNum);

% Trials
block = 1; % this is only to still display a break message if the last trial in a block is skipped due to fixation break
eyeSkip = zeros(size(trials,1),1); % trials skipped due to an eye movement, same size as trials matrix
lastFewAcc = [];
stairIdx = numel(p.stairs); % start easy
timing.startTime = GetSecs;
trialCounter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option to load a previous run from a saved workspace file (the TEMP.mat file)
% note this will overwrite most of the settings generated above
if ~isempty(workspaceFile)
    eyeFile0 = eyeFile;
    load(workspaceFile)
    fprintf('\nLOADING FROM WORKSPACE FILE:\n%s\n\n', workspaceFile)
    subjectID(end+1) = 'W';
    eyeFile = eyeFile0; % use a new eye file
%     p.forwardMask = 0; % ad only
%     block = floor(iTrial/p.nTrialsPerBlock); % ad only
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
devNum = -1;

% Present trials
while trialCounter <= nTrials
    WaitSecs(p.iti);
    
    % Initialize for eye tracking trial breaks
    if trialCounter>1
        eyeSkip(trialIdx) = stopThisTrial; % this is for the previous trial
    end
    stopThisTrial = 0;
    
    % Current trial number
    iTrial = trialCounter; % the trial we're on now
    trialIdx = trialOrder(iTrial); % the current index into trials
    trialCounter = trialCounter+1; % update the trial counter so that we will move onto the next trial, even if there is a fixation break
    
    % Get conditions for this trial
    tcCond = trials(trialIdx,targetContrastIdx);
    ts1Cond = trials(trialIdx,target1StateIdx);
    ts2Cond = trials(trialIdx,target2StateIdx);
    respInterval = p.respInterval(trials(trialIdx,respIntervalIdx));
    cueValidity = p.cueValidity(trials(trialIdx,cueValidityIdx));
    
    % Determine cued interval
    switch cueValidity
        case 1 % valid cue
            cuedInterval = respInterval;
        case -1
            cuedInterval = 3-respInterval; % 2 if 1, and 1 if 2
        case 0
            cuedInterval = 0;
    end

    % Determine target state (eg. orientation)
    switch respInterval
        case 1
            respTargetState = p.targetStates(ts1Cond);
            respTargetCond = ts1Cond;
        case 2
            respTargetState = p.targetStates(ts2Cond);
            respTargetCond = ts2Cond;
    end
    
    % Get rotations and phases
    rot = targetRotations(trialIdx,:);
    ph = targetPhases(trialIdx,:);
    
    % Select tones
    if cuedInterval==0
        cueTone = sum(p.cueTones,1)./2; % play both tones
    else
        cueTone = p.cueTones(cuedInterval,:);
    end
    respTone = p.cueTones(respInterval,:);
    
    % Select target textures
    tex1 = targetTex(tcCond,ts1Cond,ph(1));
    tex2 = targetTex(tcCond,ts2Cond,ph(2));
    
    % Select mask texture randomly
    maskIdx = randi(numel(maskTexs));
    maskTex = maskTexs(maskIdx);
    
    % Set rotation based on staircase (bypass previous)
    if p.staircase
        rotDirs = [-1 1];
        extraRot(1) = rotDirs(ts1Cond)*p.stairs(stairIdx);
        extraRot(2) = rotDirs(ts2Cond)*p.stairs(stairIdx);
        rot = rot + extraRot; % + [45 45];
        targetRotations(trialIdx,:) = rot;
    end
    
    % Store presented trial info
    trialsPresented.trials(iTrial,:) = nan(1,length(trials_headers));
    trialsPresented.trials(iTrial,1:nBaseCond) = trials(trialIdx,1:nBaseCond);
    trialsPresented.vals(iTrial).trialIdx = trialIdx;
    trialsPresented.vals(iTrial).tcCond = tcCond;
    trialsPresented.vals(iTrial).ts1Cond = ts1Cond;
    trialsPresented.vals(iTrial).ts2Cond = ts2Cond;
    trialsPresented.vals(iTrial).respInterval = respInterval;
    trialsPresented.vals(iTrial).cueValidity = cueValidity;
    trialsPresented.vals(iTrial).cuedInterval = cuedInterval;
    trialsPresented.vals(iTrial).respTargetState = respTargetState;
    trialsPresented.vals(iTrial).respTargetCond = respTargetCond;
    trialsPresented.vals(iTrial).rot = rot;
    trialsPresented.vals(iTrial).ph = ph;
    
    % Present fixation
    DrawFormattedText(window, 'x', 'center', 'center', white);
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    timeFix = Screen('Flip', window);
    
    % Check fixation hold
    if p.eyeTracking
        driftCorrected = rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});
        
        if driftCorrected
            % restart trial
            DrawFormattedText(window, 'x', 'center', 'center', white);
            drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
            timeFix = Screen('Flip', window);
        end
    end

    % Present cue
    PsychPortAudio('FillBuffer', pahandle, cueTone);
    timeCue = PsychPortAudio('Start', pahandle, [], timeFix + p.preCueDur, 1);
    % Note: waitForStart = 1 in order to return a timestamp of playback
    % More useful sound commands:
	% status = PsychPortAudio('GetStatus', pahandle);
    % soundsc(cueTone, p.Fs)
    if p.eyeTracking
        Eyelink('Message', 'EVENT_CUE');
    end
    
    % Check for eye movements
    if p.eyeTracking
        while GetSecs < timeCue + p.soas(1) - p.eyeSlack && ~stopThisTrial
            WaitSecs(.01);
            %             fixation = mod(iTrial,10); %%% for testing
            fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            [stopThisTrial trialOrder, nTrials] = fixationBreakTasks(...
                fixation, window, white*p.backgroundColor, trialOrder, iTrial, nTrials);
        end
        fixCue(iTrial) = fixation;
        if stopThisTrial
            continue
        end
    end
    
    % Present images
    if p.forwardMask(1) && ~strcmp(p.maskType,'none')
        % mask 1 - forward mask
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        Screen('DrawTexture', window, maskTex, [], imRect);
        timeMask1f = Screen('Flip', window, timeCue + p.soas(1) - p.forwardMaskSOA - slack);
        
        % blank
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeMaskBlank1f = Screen('Flip', window, timeMask1f + p.maskDur - slack);
    else
        timeMask1f = NaN;
        timeMaskBlank1f = NaN;
    end
    
    % target 1
    DrawFormattedText(window, 'x', 'center', 'center', white);
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    Screen('DrawTexture', window, tex1, [], imRect, rot(1));
    timeIm1 = Screen('Flip', window, timeCue + p.soas(1) - slack);
    if p.eyeTracking
        Eyelink('Message', 'EVENT_T1');
    end
    
    % blank
    if p.maskSOA > p.targetDur
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeBlank1 = Screen('Flip', window, timeIm1 + p.targetDur - slack);
    else
        timeBlank1 = NaN;
    end

    % mask 1 - backward mask
    if ~strcmp(p.maskType,'none') && p.backwardMask(1)
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        Screen('DrawTexture', window, maskTex, [], imRect);
        timeMask1 = Screen('Flip', window, timeIm1 + p.maskSOA - slack);
        
        % blank
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeMaskBlank1 = Screen('Flip', window, timeMask1 + p.maskDur - slack);
    else
        timeMask1 = NaN;
        timeMaskBlank1 = NaN;
    end
    
    % Check for eye movements
    if p.eyeTracking
        while GetSecs < timeCue + p.soas(2) - p.eyeSlack && ~stopThisTrial
            WaitSecs(.01);
%             fixation = mod(iTrial,10); %%% for testing
            fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            [stopThisTrial trialOrder, nTrials] = fixationBreakTasks(...
                fixation, window, white*p.backgroundColor, trialOrder, iTrial, nTrials);
        end
        fixT1(iTrial) = fixation;
        if stopThisTrial
            continue
        end
    end
    
    % mask 2 - forward mask
    % to display the T2 forward mask, the mask must be scheduled to start
    % *after* the end of the T1 backward mask
    if p.forwardMask(2) && (timeCue + p.soas(2) - p.forwardMaskSOA) > (timeMask1 + p.maskDur) && ~strcmp(p.maskType,'none')
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        Screen('DrawTexture', window, maskTex, [], imRect);
        timeMask2f = Screen('Flip', window, timeCue + p.soas(2) - p.forwardMaskSOA - slack);
        
        % blank
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeMaskBlank2f = Screen('Flip', window, timeMask2f + p.maskDur - slack);
    else
        timeMask2f = NaN;
        timeMaskBlank2f = NaN;
    end
    
    % target 2
    DrawFormattedText(window, 'x', 'center', 'center', white);
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    Screen('DrawTexture', window, tex2, [], imRect, rot(2));
    timeIm2 = Screen('Flip', window, timeCue + p.soas(2) - slack);
    if p.eyeTracking
        Eyelink('Message', 'EVENT_T2');
    end
    
    % blank
    if p.maskSOA > p.targetDur
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeBlank2 = Screen('Flip', window, timeIm2 + p.targetDur - slack);
    else
        timeBlank2 = NaN;
    end
    
    % mask 2 - backward mask
    if ~strcmp(p.maskType,'none') && p.backwardMask(2)
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        Screen('DrawTexture', window, maskTex, [], imRect);
        timeMask2 = Screen('Flip', window, timeIm2 + p.maskSOA - slack);
        
        % blank
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center', white);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeMaskBlank2 = Screen('Flip', window, timeMask2 + p.maskDur - slack);
    else
        timeMask2 = NaN;
        timeMaskBlank2 = NaN;
    end
    
    % Check for eye movements
    if p.eyeTracking
        while GetSecs < timeCue + p.respCueSOA - p.eyeSlack && ~stopThisTrial
            WaitSecs(.01);
            fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            [stopThisTrial trialOrder, nTrials] = fixationBreakTasks(...
                fixation, window, white*p.backgroundColor, trialOrder, iTrial, nTrials);
        end
        fixT2(iTrial) = fixation;
        if stopThisTrial
            continue
        end
    end
    
    % Present response cue
    PsychPortAudio('FillBuffer', pahandle, respTone);
    timeRespCue = PsychPortAudio('Start', pahandle, [], timeCue + p.respCueSOA, 1);
    if p.eyeTracking
        Eyelink('Message', 'EVENT_RESPCUE');
    end
    
    % Check for eye movements
    if p.eyeTracking
        while GetSecs < timeRespCue + p.respGoSOA - p.eyeSlack && ~stopThisTrial
            WaitSecs(.01);
            fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            [stopThisTrial trialOrder, nTrials] = fixationBreakTasks(...
                fixation, window, white*p.backgroundColor, trialOrder, iTrial, nTrials);
        end
        fixRespCue(iTrial) = fixation;
        if stopThisTrial
            continue
        end
    end
    
    % Present go cue (indicating you're allowed to make a response)
    if p.respGoSOA > 0
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center', white*p.goCueColor);
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeGoCue = Screen('Flip', window, timeRespCue + p.respGoSOA - slack);
    else
        timeGoCue = timeRespCue;
    end
    
    % Collect response
    responseKey = [];
    while isempty(responseKey) % record wrong key as missed trial
        [secs, keyCode] = KbWait(devNum);
        rt = secs - timeRespCue;
        responseKey = find(p.keyCodes==find(keyCode));
        if numel(responseKey)>1 % if more than one key was pressed simultaneously
            responseKey = [];
        end
    end
    response = p.targetStates(responseKey);
    if p.eyeTracking
        Eyelink('Message', 'TRIAL_END');
    end
    
    % Feedback
%     if response==respTargetState;
    if responseKey==respTargetCond;
        correct = 1;
        feedbackText = '+';
        feedbackColor = [0 1 0]*white;
    else
        correct = 0;
        feedbackText = '-';
        feedbackColor = [1 0 0]*white;
    end
    
    DrawFormattedText(window, feedbackText, 'center', 'center', feedbackColor);
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    timeFeedback = Screen('Flip', window);
    
    % Adjust staircase level
    if p.staircase
        [stairIdx lastFewAcc] = updateStaircase(...
            p.stairs, stairIdx, lastFewAcc, correct);
        stairValues(iTrial) = stairIdx;
        
        % Show the current threshold estimate
        reversalValues = getReversalValues(stairValues);
        if numel(reversalValues)>=5
            % average the last 5 reversals to get threshold
            threshold = mean(p.stairs(reversalValues(end-4:end)));
            fprintf('Threshold estimate = %f\n', threshold)
        else
            threshold = [];
        end
    end
   
    % Store trial info
    trials(trialIdx,cuedIntervalIdx) = cuedInterval;
    trials(trialIdx,respTargetOrientIdx) = respTargetState;
    trials(trialIdx,rtIdx) = rt;
    trials(trialIdx,responseKeyIdx) = responseKey;
    trials(trialIdx,responseIdx) = response;
    trials(trialIdx,correctIdx) = correct;
        
    % Store timing
    timing.timeFix(trialIdx,1) = timeFix;
    timing.timeCue(trialIdx,1) = timeCue;
    timing.timeMask1f(trialIdx,1) = timeMask1f;
    timing.timeMaskBlank1f(trialIdx,1) = timeMaskBlank1f;
    timing.timeIm1(trialIdx,1) = timeIm1;
    timing.timeBlank1(trialIdx,1) = timeBlank1;
    timing.timeMask1(trialIdx,1) = timeMask1;
    timing.timeMaskBlank1(trialIdx,1) = timeMaskBlank1;
    timing.timeMask2f(trialIdx,1) = timeMask2f;
    timing.timeMaskBlank2f(trialIdx,1) = timeMaskBlank2f;
    timing.timeIm2(trialIdx,1) = timeIm2;
    timing.timeBlank2(trialIdx,1) = timeBlank2;
    timing.timeMask2(trialIdx,1) = timeMask2;
    timing.timeMaskBlank2(trialIdx,1) = timeMaskBlank2;
    timing.timeRespCue(trialIdx,1) = timeRespCue;
    timing.timeGoCue(trialIdx,1) = timeGoCue;
    timing.timeFeedback(trialIdx,1) = timeFeedback;

    % Store presented trial info
    trialsPresented.trials(iTrial,nBaseCond+1:end) = trials(trialIdx,nBaseCond+1:end);

    
    save('data/TEMP') % saves the workspace on each trial
    
    if mod(iTrial,p.nTrialsPerBlock)==0 || iTrial==nTrials || iTrial>p.nTrialsPerBlock*block        
        % Save workspace if practice block (since we might quit after just
        % one block)
        if ~isempty(strfind(subjectID,'PRAC'))
            save(sprintf('data/WORKSPACE_%s', subjectID))
        end
        
        % Calculate block accuracy
        blockStartTrial = (iTrial/p.nTrialsPerBlock)*p.nTrialsPerBlock - p.nTrialsPerBlock + 1;
        if blockStartTrial < 0 % we are doing less than one block
            blockStartTrial = 1;
        end
        inBlock = zeros(size(trials,1),1);
        inBlock(trialOrder(blockStartTrial:iTrial)) = 1;
        completedInBlock = inBlock & ~eyeSkip;
        nTrialsCompletedInBlock = nnz(completedInBlock)
        blockAcc = mean(trials(completedInBlock,correctIdx))
        
        accMessage = sprintf('Accuracy: %d%%', round(blockAcc*100));
        blockMessage = sprintf('%s You''ve completed %d of %d blocks.', highpraise, block, ceil(nTrials/p.nTrialsPerBlock));
        if iTrial==nTrials
            keyMessage = '';
        else
            keyMessage = 'Press any key to go on.';
        end
        
        breakMessage = sprintf('%s\n%s\n\n%s', blockMessage, accMessage, keyMessage);
        DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window);
        WaitSecs(1);
        if iTrial < nTrials
            KbWait(devNum);
        end
        
        block = block+1; % keep track of block for block message only
    end
end
timing.endTime = GetSecs;

DrawFormattedText(window, 'All done! Thanks for your effort!', 'center', 'center', [1 1 1]*white);
Screen('Flip', window);
WaitSecs(2);

%% Store expt info
expt.subjectID = subjectID;
expt.p = p;
expt.timing = timing;
expt.trialOrder = trialOrder;
expt.trials_headers = trials_headers;
expt.trials = trials;
expt.targetRotations = targetRotations;
expt.targetPhases = targetPhases;
expt.trialsPresented = trialsPresented;

if p.staircase
    expt.staircase.stairValues = stairValues;
    expt.staircase.reversalValues = reversalValues;
    expt.staircase.threshold = threshold;
end

if p.eyeTracking
    expt.eye.fixCue = fixCue;
    expt.eye.fixT1 = fixT1;
    expt.eye.fixT2 = fixT2;
    expt.eye.fixRespCue = fixRespCue;
end

%% Calculate more timing things
expt.timing.dur.im1 = expt.timing.timeBlank1 - expt.timing.timeIm1;
expt.timing.dur.im2 = expt.timing.timeBlank2 - expt.timing.timeIm2;
expt.timing.dur.cueIm1SOA = expt.timing.timeIm1 - expt.timing.timeCue;
expt.timing.dur.cueIm2SOA = expt.timing.timeIm2 - expt.timing.timeCue;
expt.timing.dur.im1Im2SOA = expt.timing.timeIm2 - expt.timing.timeIm1;

%% Analyze and save data
[expt results] = rd_analyzeTemporalAttention(expt, saveData, saveFigs, plotTimingFigs, saveTimingFigs);

%% Save eye data and shut down the eye tracker
if p.eyeTracking
    rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});
    
    % rename eye file
    eyeFileFull = sprintf('%s/%s_TemporalAttention_%s.edf', eyeDataDir, subjectID, datestr(now, 'yyyymmdd'));
    copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)
end

%% Clean up
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll')
