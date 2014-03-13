function rd_temporalAttention(subjectID)

%% Setup
if nargin==0
    subjectID = 'test';
end

saveData = 1;

p = temporalAttentionParams;

slack = p.refRate/2;

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

%% Keyboard
% Find keyboard device number
devNum = findKeyboardDevNumsAtLocationNYU(p.testingLocation);
if isempty(devNum)
    error('Could not find Keypad!')
end

%% Sound
% Perform basic initialization of the sound driver
InitializePsychSound(1); % 1 for precise timing

% Open audio device for low-latency output
reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel

%% Screen
% Set up window and textures
screenNumber = max(Screen('Screens'));
% [window rect] = Screen('OpenWindow', screenNumber);
[window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]);
white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
[cx cy] = RectCenter(rect);
Screen('TextSize', window, 24);
Screen('TextColor', window, [1 1 1]*white);
Screen('TextFont', window, p.font);

% Check screen size
[sw, sh] = Screen('WindowSize', window); % height and width of screen (px)
if ~all([sw sh] == p.screenRes) && ~strcmp(subjectID,'test')
    error('Screen resolution is different from requested!')
end

% Check refresh rate
flipInterval = Screen('GetFlipInterval', window); % frame duration (s)
if abs(flipInterval - p.refRate) > 0.001
    error('Refresh rate is different from requested!')
end

% Check font
if ~strcmp(p.font, Screen('TextFont', window))
    error('Font was not set to requested: %s', p.font)
end

% Load calibration file
switch p.testingLocation
    case 'CarrascoL1'
        load('../Displays/0001_james_TrinitonG520_1280x960_57cm_Input1_140129.mat');
        Screen('LoadNormalizedGammaTable', window, repmat(calib.table,1,3));
    otherwise
        fprintf('\nNot loading gamma table ...\n')
end

%% Make stimuli
% Calculate stimulus dimensions (px) and position
imPos = round(ang2pix(p.imPos, p.screenSize(1), p.screenRes(1), p.viewDist, 'radial')); % from screen center
targetSize = round(ang2pix(p.targetSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
pixelsPerDegree = round(ang2pix(1, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));

% Make target gratings
for iC = 1:numel(p.targetContrasts)
    targetContrast = p.targetContrasts(iC);
    for iTO = 1:numel(p.targetOrientations)
        targetOrientation = p.targetOrientations(iTO);
        
        % make big grating
        t = buildColorGrating(pixelsPerDegree, p.imSize, ...
            p.spatialFrequency, targetOrientation, 0, targetContrast, 0, 'bw');
        
        % mask with annulus
        target{iC,iTO} = maskWithGaussian(t, size(t,1), targetSize);
    end
end

% Make mask
switch p.maskType
    case 'whitenoise'
        mask = (rand(size(t))-0.5)*p.maskContrast + 0.5;
    case 'verticalgrating'
        m = buildColorGrating(pixelsPerDegree, p.imSize, ...
            p.spatialFrequency, 0, 0, p.maskContrast, 0, 'bw');
        mask = maskWithGaussian(m, size(m,1), targetSize);
    case 'crossedgratings'
        for i=1:numel(p.targetOrientations)
            m(:,:,i) = buildColorGrating(pixelsPerDegree, p.imSize, ...
                p.spatialFrequency, p.targetOrientations(i), 0, p.maskContrast, 0, 'bw');
        end
        m = sum(m,3)./numel(p.targetOrientations);
        mask = maskWithGaussian(m, size(m,1), targetSize);
    case 'filterednoise'
        m = rand(size(t));
        f = lowpassfilter(size(m), 0.15, 5);
        m = imifft(imfft(m).*f);
%         m = (m-min(m(:))) ./ (max(m(:)-min(m(:)))); % rescale
        m = histeq(m); % basically increases the contrast
        mask = (m-0.5)*p.maskContrast + 0.5;
    otherwise
        error('maskType not recognized')
end

% Make textures
for iC = 1:numel(p.targetContrasts)
    for iTO = 1:numel(p.targetOrientations)
        targetTex(iC,iTO) = Screen('MakeTexture', window, target{iC,iTO}*white);
    end
end

maskTex = Screen('MakeTexture', window, mask*white);

% Make the rects for placing the images
imSize = size(target{1,1},1);
imRect = CenterRectOnPoint([0 0 imSize imSize], cx+imPos(1), cy+imPos(2));

%% Generate trials
% Construct trials matrix
trials_headers = {'targetContrast','respInterval','cueValidity',...
    'target1Orient','target2Orient','cuedInterval','respTargetOrient',...
    'rt','responseKey','response','correct'};

% make sure column indices match trials headers
targetContrastIdx = strcmp(trials_headers,'targetContrast');
respIntervalIdx = strcmp(trials_headers,'respInterval');
cueValidityIdx = strcmp(trials_headers,'cueValidity');
target1OrientIdx = strcmp(trials_headers,'target1Orient');
target2OrientIdx = strcmp(trials_headers,'target2Orient');
cuedIntervalIdx = strcmp(trials_headers,'cuedInterval');
respTargetOrientIdx = strcmp(trials_headers,'respTargetOrient');
rtIdx = strcmp(trials_headers,'rt');
responseKeyIdx = strcmp(trials_headers,'responseKey');
responseIdx = strcmp(trials_headers,'response');
correctIdx = strcmp(trials_headers,'correct');

% full factorial design
trials = fullfact([numel(p.targetContrasts) ...
    numel(p.respInterval) ...
    numel(p.cueValidityFactor) ...
    numel(p.targetOrientations) ...
    numel(p.targetOrientations)]);

% set cue validity condition according to the cueValidityFactor (potentially unequal proportion of trials in each condition)
cueValidityTrials = trials(:,cueValidityIdx);
for i = 1:numel(p.cueValidityFactor)   
    trials(cueValidityTrials==i,cueValidityIdx) = p.cueValidityFactor(i);
end

% repeat trials matrix according to nReps of all conditions
trials = repmat(trials, p.nReps, 1);
nTrials = size(trials,1);

% Choose order of trial presentation
trialOrder = randperm(nTrials);

%% Present trials
% Show fixation and wait for a button press
Screen('FillRect', window, white*p.backgroundColor);
DrawFormattedText(window, 'x', 'center', 'center');
Screen('Flip', window);
KbWait(devNum);

% Trials
timing.startTime = GetSecs;
for iTrial = 1:nTrials
    WaitSecs(p.iti);
    
    % Current trial number
    trialIdx = trialOrder(iTrial);
    
    % Get conditions for this trial
    tcCond = trials(trialIdx,targetContrastIdx);
    to1Cond = trials(trialIdx,target1OrientIdx);
    to2Cond = trials(trialIdx,target2OrientIdx);
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

    % Determine target orientation
    switch respInterval
        case 1
            respTargetOrientation = p.targetOrientations(to1Cond);
        case 2
            respTargetOrientation = p.targetOrientations(to2Cond);
    end
    
    % Select tones and textures
    if cuedInterval==0
        cueTone = sum(p.cueTones,1); % play both tones
    else
        cueTone = p.cueTones(cuedInterval,:);
    end
    respTone = p.cueTones(respInterval,:);
    
    tex1 = targetTex(tcCond,to1Cond);
    tex2 = targetTex(tcCond,to2Cond);
    
    % Present fixation
    DrawFormattedText(window, 'x', 'center', 'center', [1 1 1]*white);
    timeFix = Screen('Flip', window);
    
    % Present cue
    PsychPortAudio('FillBuffer', pahandle, cueTone);
    timeCue = PsychPortAudio('Start', pahandle, [], timeFix + p.preCueDur, 1);
    % Note: waitForStart = 1 in order to return a timestamp of playback
    % More useful sound commands:
	% status = PsychPortAudio('GetStatus', pahandle);
    % soundsc(cueTone, p.Fs)
    
    % Present images
    % target 1
    Screen('DrawTexture', window, tex1, [], imRect);
    DrawFormattedText(window, 'x', 'center', 'center');
    timeIm1 = Screen('Flip', window, timeCue + p.soas(1) - slack);
    
    % blank
    Screen('FillRect', window, white*p.backgroundColor);
    DrawFormattedText(window, 'x', 'center', 'center');
    timeBlank1 = Screen('Flip', window, timeIm1 + p.targetDur - slack);
    
    % mask 1
    Screen('DrawTexture', window, maskTex, [], imRect);
    DrawFormattedText(window, 'x', 'center', 'center');
    timeMask1 = Screen('Flip', window, timeIm1 + p.maskSOA - slack);

    % blank
    Screen('FillRect', window, white*p.backgroundColor);
    DrawFormattedText(window, 'x', 'center', 'center');
    timeMaskBlank1 = Screen('Flip', window, timeMask1 + p.maskDur - slack);
    
    % target 2
    Screen('DrawTexture', window, tex2, [], imRect, 0);
    DrawFormattedText(window, 'x', 'center', 'center');
    timeIm2 = Screen('Flip', window, timeCue + p.soas(2) - slack);
    
    % blank
    Screen('FillRect', window, white*p.backgroundColor);
    DrawFormattedText(window, 'x', 'center', 'center');
    timeBlank2 = Screen('Flip', window, timeIm2 + p.targetDur - slack);
    
    % mask 2
    Screen('DrawTexture', window, maskTex, [], imRect);
    DrawFormattedText(window, 'x', 'center', 'center');
    timeMask2 = Screen('Flip', window, timeIm2 + p.maskSOA - slack);

    % blank
    Screen('FillRect', window, white*p.backgroundColor);
    DrawFormattedText(window, 'x', 'center', 'center');
    timeMaskBlank2 = Screen('Flip', window, timeMask2 + p.maskDur - slack);
    
    % Present response cue
    PsychPortAudio('FillBuffer', pahandle, respTone);
    timeRespCue = PsychPortAudio('Start', pahandle, [], timeCue + p.respCueSOA, 1);
    
    % Collect response
    responseKey = [];
    while isempty(responseKey) % record wrong key as missed trial
        [secs, keyCode] = KbWait(devNum);
        rt = secs - timeRespCue;
        responseKey = find(p.keyCodes==find(keyCode));
    end
    response = p.targetOrientations(responseKey);
    
    % Feedback
    if response==respTargetOrientation;
        correct = 1;
        feedbackText = '+';
        feedbackColor = [0 1 0]*white;
    else
        correct = 0;
        feedbackText = '-';
        feedbackColor = [1 0 0]*white;
    end
    
    DrawFormattedText(window, feedbackText, 'center', 'center', feedbackColor);
    timeFeedback = Screen('Flip', window);
   
    % Store trial info
    trials(trialIdx,cuedIntervalIdx) = cuedInterval;
    trials(trialIdx,respTargetOrientIdx) = respTargetOrientation;
    trials(trialIdx,rtIdx) = rt;
    trials(trialIdx,responseKeyIdx) = responseKey;
    trials(trialIdx,responseIdx) = response;
    trials(trialIdx,correctIdx) = correct;
    
    % Store timing
    timing.timeFix(iTrial,1) = timeFix;
    timing.timeCue(iTrial,1) = timeCue;
    timing.timeIm1(iTrial,1) = timeIm1;
    timing.timeBlank1(iTrial,1) = timeBlank1;
    timing.timeIm2(iTrial,1) = timeIm2;
    timing.timeBlank2(iTrial,1) = timeBlank2;
    timing.timeRespCue(iTrial,1) = timeRespCue;
    timing.timeFeedback(iTrial,1) = timeFeedback;
    
    save('data/TEMP') % saves the workspace on each trial
end
timing.endTime = GetSecs;

% Show end of block feedback
% acc = mean(trials(:,correctIdx));

% Store expt info
expt.p = p;
expt.timing = timing;
expt.trialOrder = trialOrder;
expt.trials_headers = trials_headers;
expt.trials = trials;

% Clean up
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll')

% Analyze data
for iCV = 1:numel(p.cueValidity)
    w = trials(:,cueValidityIdx)==iCV;
    totals.all{iCV} = trials(w,:);
    
    totals.means(iCV,:) = mean(totals.all{iCV},1);
    totals.stds(iCV,:) = std(totals.all{iCV},0,1);
    totals.stes(iCV,:) = totals.stds(iCV,:)./sqrt(size(totals.all{iCV},1));
end

accMean = totals.means(:,correctIdx);
accSte = totals.stes(:,correctIdx);

rtMean = totals.means(:,rtIdx);
rtSte = totals.stes(:,rtIdx);

% Plot figs
figure
errorbar(p.cueValidity, accMean, accSte, '.k')
xlabel('cue validity')
ylabel('acc')

figure
errorbar(p.cueValidity, rtMean, rtSte, '.k')
xlabel('cue validity')
ylabel('rt')

% Store data
results.totals = totals;
results.accMean = accMean;
results.accSte = accSte;
results.whenSaved = datestr(now);

% Save data
if saveData
    fileName = sprintf('data/%s_temporalAttention_%s.mat', subjectID, datestr(now, 'yyyymmdd'));
    save(fileName, 'expt', 'results')
end


