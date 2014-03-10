function rd_temporalAttention(subjectID)

%% Setup
if nargin==0
    subjectID = 'test';
end

saveData = 1;

p = temporalAttentionParams;

slack = p.refRate/2;

%% Keyboard and screen
% Find keyboard device number
devNum = findKeyboardDevNumsAtLocationNYU(p.testingLocation);
if isempty(devNum)
    error('Could not find Keypad!')
end

% Set up window and textures
screenNumber = max(Screen('Screens'));
% [win rect] = Screen('OpenWindow', screenNumber);
[win rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]);
white = WhiteIndex(win);  % Retrieves the CLUT color code for white.
[cx cy] = RectCenter(rect);

% Check screen size
[sw, sh] = Screen('WindowSize', win); % height and width of screen (px)
if ~all([sw sh] == p.screenRes) && ~strcmp(subjectID,'test')
    error('Screen resolution is different from requested!')
end

% Check refresh rate
flipInterval = Screen('GetFlipInterval', win); % frame duration (s)
if abs(flipInterval - p.refRate) > 0.001
    error('Refresh rate is different from requested!')
end

% Load calibration file
switch p.testingLocation
    case 'CarrascoL1'
        load('../Displays/0001_james_TrinitonG520_1280x960_57cm_Input1_140129.mat');
        Screen('LoadNormalizedGammaTable', win, repmat(calib.table,1,3));
    otherwise
        fprintf('\nNot loading gamma table ...\n')
end

%% Make stimuli
% Calculate stimulus dimensions (px) and position
imPos = round(ang2pix(p.imPos, p.screenSize(1), p.screenRes(1), p.viewDist, 'radial')); % from screen center
targetSize = round(ang2pix(p.targetSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
fixSize = round(ang2pix(p.fixSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
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

% Make textures
for iC = 1:numel(p.targetContrasts)
    for iTO = 1:numel(p.targetOrientations)
        targetTex(iC,iTO) = Screen('MakeTexture', win, target{iC,iTO}*white);
    end
end

% Make the rects for placing the images
imSize = size(target{1,1},1);
imRect = CenterRectOnPoint([0 0 imSize imSize], cx+imPos(1), cy+imPos(2));
fixRect = CenterRectOnPoint([0 0 fixSize fixSize], cx, cy);

%% Generate trials
% Construct trials matrix
trials_headers = {'targetContrast','cuedInterval','cueValidity','target1Orient','target2Orient','rt','responseKey','response','correct'};

targetContrastIdx = find(strcmp(trials_headers,'targetContrast'));
cuedIntervalIdx = find(strcmp(trials_headers,'cuedInterval'));
cueValidityIdx = find(strcmp(trials_headers,'cueValidity'));
target1OrientIdx = find(strcmp(trials_headers,'target1Orient'));
target2OrientIdx = find(strcmp(trials_headers,'target2Orient'));

% full factorial design
trials = fullfact([numel(p.targetContrasts) ...
    numel(p.cuedInterval) ...
    numel(p.cueValidityFactor) ...
    numel(p.targetOrientations) ...
    numel(p.targetOrientations)]);

% set cue validity to be just 1 or 2 (in potentially unequal proportions)
cueValidityTrials = trials(:,cueValidityIdx);
for i = 1:numel(p.cueValidityFactor)   
    trials(cueValidityTrials==i,cueValidityIdx) = p.cueValidityFactor(i);
end

trials = repmat(trials, p.nReps, 1);
nTrials = size(trials,1);

% Choose order of trial presentation
trialOrder = randperm(nTrials);

%% Present trials
% Show fixation and wait for a button press
Screen('FillRect', win, white*p.backgroundColor);
Screen('FillRect', win, [0 0 0], fixRect);
vbl = Screen('Flip', win);
KbWait(devNum);

% Trials
timing.startTime = GetSecs;
for iTrial = 1:nTrials
    % Current trial number
    trialIdx = trialOrder(iTrial);
    
    % Get conditions for this trial
    tcCond = trials(trialIdx,targetContrastIdx);
    to1Cond = trials(trialIdx,target1OrientIdx);
    to2Cond = trials(trialIdx,target2OrientIdx);
    targetContrast = p.targetContrasts(tcCond);
    cuedInterval = p.cuedInterval(trials(trialIdx,cuedIntervalIdx));
    cueValidity = p.cuedInterval(trials(trialIdx,cueValidityIdx));
    
    % Determine response interval and target orientation
    if cueValidity==1 % valid cue
        respInterval = cuedInterval;
    else
        respInterval = 3-cuedInterval; % 2 if 1, and 1 if 2
    end

    if respInterval==1
        respTargetOrientation = p.targetOrientations(to1Cond);
    else
        respTargetOrientation = p.targetOrientations(to2Cond);
    end
    
    % Select tones and textures
    cueTone = p.cueTones(cuedInterval,:);
    respTone = p.cueTones(respInterval,:);
    
    tex1 = targetTex(tcCond,to1Cond);
    tex2 = targetTex(tcCond,to2Cond);
    
    % Present cue
    %%% insert tone cue here %%%
    soundsc(cueTone)
    Screen('FillRect', win, [255 255 255], fixRect);
    timeCue = Screen('Flip', win);
    
    % Present images
    Screen('DrawTexture', win, tex1, [], imRect);
    Screen('FillRect', win, [255 255 255], fixRect);
    timeIm1 = Screen('Flip', win, timeCue + p.soas(1) - slack);
    
    Screen('FillRect', win, white*p.backgroundColor);
    Screen('FillRect', win, [255 255 255], fixRect);
    timeBlank1 = Screen('Flip', win, timeIm1 + p.targetDur - slack);
    
    Screen('DrawTexture', win, tex2, [], imRect);
    Screen('FillRect', win, [255 255 255], fixRect);
    timeIm2 = Screen('Flip', win, timeCue + p.soas(2) - slack);
    
    Screen('FillRect', win, white*p.backgroundColor);
    Screen('FillRect', win, [255 255 255], fixRect);
    timeBlank2 = Screen('Flip', win, timeIm2 + p.targetDur - slack);
    
    %%% insert tone respone cue here %%%
    soundsc(respTone)
    Screen('FillRect', win, [0 0 255], fixRect);
    timeRespCue = Screen('Flip', win, timeCue + p.respCueSOA - slack);
    
    % Collect response
    [secs, keyCode] = KbWait(devNum);
    rt = secs - timeRespCue;
    
    % Feedback
    responseKey = find(p.keyCodes==find(keyCode));
    response = p.targetOrientations(responseKey);
    if response==respTargetOrientation;
        correct = 1;
        feedbackText = '+';
        feedbackColor = [0 1 0]*white;
    else
        correct = 0;
        feedbackText = '-';
        feedbackColor = [1 0 0]*white;
    end
    
    DrawFormattedText(win, feedbackText, cx, cy, feedbackColor)
    timeFeedback = Screen('Flip', win);
   
    % Store trial info
    trials(trialIdx,6) = respInterval;
    trials(trialIdx,7) = respTargetOrientation;
    trials(trialIdx,8) = rt;
    trials(trialIdx,9) = responseKey;
    trials(trialIdx,10) = response;
    trials(trialIdx,11) = correct;
    
    % Store timing
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
% acc = mean(trials(:,6));

% Store expt info
expt.p = p;
expt.timing = timing;
expt.trialOrder = trialOrder;
expt.trials_headers = trials_headers;
expt.trials = trials;

% Clean up
Screen('CloseAll')

% Analyze data
for iSOA = 1:numel(p.soas)
    w = trials(:,1)==iSOA;
    totals.all(:,:,iSOA) = trials(w,:);
end

totals.means = squeeze(mean(totals.all,1))';
totals.stds = squeeze(std(totals.all,0,1))';
totals.stes = totals.stds./sqrt(size(totals.all,1));

accMean = totals.means(:,8);
accSte = totals.stes(:,8);

% Plot figs
figure
errorbar(p.soas, accMean, accSte, '.-k')

% figure
% scatter(trials(:,3), trials(:,4))

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


