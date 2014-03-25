function rd_temporalAttention(subjectID)

%% Setup
if nargin==0
    subjectID = 'test';
end

saveData = 1;
saveFigs = 1;

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
[window rect] = Screen('OpenWindow', screenNumber);
% [window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]);
white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
[cx cy] = RectCenter(rect);
Screen('TextSize', window, 24);
Screen('TextColor', window, white);
Screen('TextFont', window, p.font);

% Set resolution and refresh rate
if strcmp(p.testingLocation, 'CarrascoL1')
    SetResolution(screenNumber, p.screenRes(1), p.screenRes(2), round(1/p.refRate));
end

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
                
                % make big grating
                t = buildColorGrating(pixelsPerDegree, p.imSize, ...
                    spatialFrequency, targetOrientation, 0, targetContrast, 0, 'bw');
        end

        % mask with annulus
        target{iC,iT} = maskWithGaussian(t, size(t,1), targetSize);
    end
end

% Make mask
switch p.maskType
    case 'none'
        mask = ones(size(t))*p.backgroundColor;
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

%% Make textures
for iC = 1:numel(p.targetContrasts)
    for iT = 1:numel(p.targetStates)
        targetTex(iC,iT) = Screen('MakeTexture', window, target{iC,iT}*white);
    end
end

maskTex = Screen('MakeTexture', window, mask*white);

% Make the rects for placing the images
imSize = size(target{1,1},1);
imRect = CenterRectOnPoint([0 0 imSize imSize], cx+imPos(1), cy+imPos(2));
phRect = imRect + [-1 -1 1 1]*p.phLineWidth;

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

% generate target rotations
% these are applied over and above any initial orientation
nTrials0 = size(trials,1);
switch p.rotateTarget
    case 'none'
        targetRotations = zeros(nTrials0,2);
    case 'random'
        targetRotations = 360*rand(nTrials0,2);
    case 'rotT2'
        targetRotations = [zeros(nTrials0,1) 90*ones(nTrials0,1)];
    case 'cb'
        % if we want to counterbalance rotations across trials, we need to
        % have 4x the original number of trials
        targetRotations(:,1) = repmat([zeros(nTrials0,1); 90*ones(nTrials0,1)],2,1);
        targetRotations(:,2) = [zeros(nTrials0*2,1); 90*ones(nTrials0*2,1)];
        trials = repmat(trials,4,1);
    otherwise
        error('p.rotateTarget not recognized')
end
        
% set cue validity condition according to the cueValidityFactor (potentially unequal proportion of trials in each condition)
cueValidityTrials = trials(:,cueValidityIdx);
for i = 1:numel(p.cueValidityFactor)   
    trials(cueValidityTrials==i,cueValidityIdx) = p.cueValidityFactor(i);
end

% repeat trials matrix according to nReps of all conditions
trials = repmat(trials, p.nReps, 1);
nTrials = size(trials,1);

% determine blocks
nBlocks = nTrials/p.nTrialsPerBlock;
fprintf('\n%1.2f blocks\n\n', nBlocks)

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
        case 2
            respTargetState = p.targetStates(ts2Cond);
    end
    
    % Select tones and textures
    if cuedInterval==0
        cueTone = sum(p.cueTones,1); % play both tones
    else
        cueTone = p.cueTones(cuedInterval,:);
    end
    respTone = p.cueTones(respInterval,:);
    
    tex1 = targetTex(tcCond,ts1Cond);
    tex2 = targetTex(tcCond,ts2Cond);

    % Get rotations
    rot = targetRotations(trialIdx,:);
    
    % Present fixation
    DrawFormattedText(window, 'x', 'center', 'center', white);
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
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
    DrawFormattedText(window, 'x', 'center', 'center');
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    Screen('DrawTexture', window, tex1, [], imRect, rot(1));
    timeIm1 = Screen('Flip', window, timeCue + p.soas(1) - slack);
    
    % blank
    if p.maskSOA > p.targetDur
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center');
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeBlank1 = Screen('Flip', window, timeIm1 + p.targetDur - slack);
    else
        timeBlank1 = NaN;
    end
    
    % mask 1
    DrawFormattedText(window, 'x', 'center', 'center');
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    Screen('DrawTexture', window, maskTex, [], imRect);
    timeMask1 = Screen('Flip', window, timeIm1 + p.maskSOA - slack);

    % blank
    Screen('FillRect', window, white*p.backgroundColor);
    DrawFormattedText(window, 'x', 'center', 'center');
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    timeMaskBlank1 = Screen('Flip', window, timeMask1 + p.maskDur - slack);
    
    % target 2
    DrawFormattedText(window, 'x', 'center', 'center');
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    Screen('DrawTexture', window, tex2, [], imRect, rot(2));
    timeIm2 = Screen('Flip', window, timeCue + p.soas(2) - slack);
    
    % blank
    if p.maskSOA > p.targetDur
        Screen('FillRect', window, white*p.backgroundColor);
        DrawFormattedText(window, 'x', 'center', 'center');
        drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
        timeBlank2 = Screen('Flip', window, timeIm2 + p.targetDur - slack);
    else
        timeBlank2 = NaN;
    end
    
    % mask 2
    DrawFormattedText(window, 'x', 'center', 'center');
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    Screen('DrawTexture', window, maskTex, [], imRect);
    timeMask2 = Screen('Flip', window, timeIm2 + p.maskSOA - slack);

    % blank
    Screen('FillRect', window, white*p.backgroundColor);
    DrawFormattedText(window, 'x', 'center', 'center');
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
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
    response = p.targetStates(responseKey);
    
    % Feedback
    if response==respTargetState;
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
   
    % Store trial info
    trials(trialIdx,cuedIntervalIdx) = cuedInterval;
    trials(trialIdx,respTargetOrientIdx) = respTargetState;
    trials(trialIdx,rtIdx) = rt;
    trials(trialIdx,responseKeyIdx) = responseKey;
    trials(trialIdx,responseIdx) = response;
    trials(trialIdx,correctIdx) = correct;
        
    % Store timing
    timing.timeFix(iTrial,1) = timeFix;
    timing.timeCue(iTrial,1) = timeCue;
    timing.timeIm1(iTrial,1) = timeIm1;
    timing.timeBlank1(iTrial,1) = timeBlank1;
    timing.timeMask1(iTrial,1) = timeMask1;
    timing.timeMaskBlank1(iTrial,1) = timeMaskBlank1;
    timing.timeIm2(iTrial,1) = timeIm2;
    timing.timeBlank2(iTrial,1) = timeBlank2;
    timing.timeMask2(iTrial,1) = timeMask2;
    timing.timeMaskBlank2(iTrial,1) = timeMaskBlank2;
    timing.timeRespCue(iTrial,1) = timeRespCue;
    timing.timeFeedback(iTrial,1) = timeFeedback;
    
    save('data/TEMP') % saves the workspace on each trial
    
    if mod(iTrial,p.nTrialsPerBlock)==0 && iTrial<nTrials
        DrawFormattedText(window, 'Break time!\n\nPress any key to go on.', 'center', 'center', [1 1 1]*white);
        Screen('Flip', window);
        WaitSecs(1);
        KbWait(devNum);
    end
end
timing.endTime = GetSecs;

% Show end of block feedback
% acc = mean(trials(:,correctIdx));

% Store expt info
expt.subjectID = subjectID;
expt.p = p;
expt.timing = timing;
expt.trialOrder = trialOrder;
expt.trials_headers = trials_headers;
expt.trials = trials;
expt.targetRotations = targetRotations;

% Clean up
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll')

% Analyze data
for iRI = 1:numel(p.respInterval)
    for iCV = 1:numel(p.cueValidity)
        for iTC = 1:numel(p.targetContrasts)
            w = trials(:,respIntervalIdx)==iRI & trials(:,cueValidityIdx)==iCV & trials(:,targetContrastIdx)==iTC;
            totals.all{iCV,iRI}(:,:,iTC) = trials(w,:);
            
            totals.means{iRI}(iCV,:,iTC) = mean(totals.all{iCV,iRI}(:,:,iTC),1);
            totals.stds{iRI}(iCV,:,iTC) = std(totals.all{iCV,iRI}(:,:,iTC),0,1);
            totals.stes{iRI}(iCV,:,iTC) = totals.stds{iRI}(iCV,:,iTC)./sqrt(size(totals.all{iCV,iRI}(:,:,iTC),1));
        end
    end
end

for iRI = 1:numel(p.respInterval)
    accMean{iRI} = squeeze(totals.means{iRI}(:,correctIdx,:)); % [validity x contrast]
    accSte{iRI} = squeeze(totals.stes{iRI}(:,correctIdx,:));
    
    rtMean{iRI} = squeeze(totals.means{iRI}(:,rtIdx,:));
    rtSte{iRI} = squeeze(totals.stes{iRI}(:,rtIdx,:));
end

% Store data
results.totals = totals;
results.accMean = accMean;
results.accSte = accSte;
results.rtMean = rtMean;
results.rtSte = rtSte;
results.whenSaved = datestr(now);

% Save data
if saveData
    fileName = sprintf('data/%s_TemporalAttention_%s.mat', subjectID, datestr(now, 'yyyymmdd'));
    save(fileName, 'expt', 'results')
end

% Plot figs
intervalNames = {'early','late'};
accLims = [0.2 1];
rtLims = [0.3 1.6];
contrastLims = [p.targetContrasts(1)-0.05 p.targetContrasts(end)+0.05];

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(contrastLims, [0.5 0.5], '--k');
    p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
        accMean{iRI}', accSte{iRI}', '.', 'MarkerSize', 20);
    xlabel('contrast')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(accLims)
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
        rtMean{iRI}', rtSte{iRI}', '.', 'MarkerSize', 20)
    xlabel('contrast')
    ylabel('rt')
    legend(num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(rtLims)
    box off
end

% Save figs
if saveFigs
    figNames = {'acc','rt'};
    rd_saveAllFigs(fig, figNames, [subjectID '_TemporalAttention'])
end

