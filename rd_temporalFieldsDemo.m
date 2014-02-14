% rd_temporalFieldsDemo.m

p = temporalFieldsParams;

slack = p.refRate/2;

% Find keyboard device number
devNum = findKeyboardDevNumsAtLocationNYU(p.testingLocation);
if isempty(devNum)
    error('Could not find Keypad!')
end

% Sound for starting trials
v = 1:1000;
envelope = [0:.05:1 ones(1,1000-42) 1:-.05:0];
startSound = (0.25 * sin(3*pi*v/30)).*envelope;

% Set up window and textures
screenNumber = max(Screen('Screens'));
[win rect] = Screen('OpenWindow', screenNumber);
% [win rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]);
white = WhiteIndex(win);  % Retrieves the CLUT color code for white.
[cx cy] = RectCenter(rect);

% Check screen size
[sw, sh] = Screen('WindowSize', win); % height and width of screen (px)
if ~all([sw sh] == p.screenRes)
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

% Calculate stimulus dimensions (px) and position
imPos = round(ang2pix(p.imPos, p.screenSize(1), p.screenRes(1), p.viewDist, 'radial')); % from screen center
targetSize = round(ang2pix(p.targetSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
surroundSize = round(ang2pix(p.surroundSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
fixSize = round(ang2pix(p.fixSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
pixelsPerDegree = round(ang2pix(1, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));

% Make target and surround gratings
% Make big gratings that will later be masked
t = buildColorGrating(pixelsPerDegree, p.imSize, ...
    p.spatialFrequency, p.orientation, 0, p.targetContrast, 0, 'bw');
s = buildColorGrating(pixelsPerDegree, p.imSize, ...
    p.spatialFrequency, p.orientation, 0, p.surroundContrast, 0, 'bw');
blank = buildColorGrating(pixelsPerDegree, p.imSize, ...
    p.spatialFrequency, p.orientation, 0, 0, 0, 'bw');

% Mask with annulus
target = maskWithGaussian(t, size(t,1), targetSize);
surround = maskWithGaussian(s, size(s,1), surroundSize);

% Make textures
targetTex = Screen('MakeTexture', win, target*white);
surroundTex = Screen('MakeTexture', win, surround*white);
blankTex = Screen('MakeTexture', win, blank*white);

% Make the rects for placing the images
imSize = size(target,1);
imRect = CenterRectOnPoint([0 0 imSize imSize], cx+imPos(1), cy+imPos(2));
fixRect = CenterRectOnPoint([0 0 fixSize fixSize], cx, cy);

% Construct trials matrix
trials_headers = {'soaCond','targetPresAbsCond','soa','actualSOA','rt','responseKey','response','correct','jitter'};

trials = fullfact([numel(p.soas) numel(p.targetPresAbs)]);
trials = repmat(trials, p.nReps, 1);
nTrials = size(trials,1);

%%% for debugging, just show all the stimuli
% trials(:,2) = 1;

% Choose order of trial presentation
trialOrder = randperm(nTrials);

% Show fixation and wait for a button press
Screen('FillRect', win, white*p.backgroundColor);
Screen('FillRect', win, [0 0 0], fixRect);
vbl = Screen('Flip', win);
KbWait(devNum);

% Trials
tic
for iTrial = 1:nTrials
    % get conditions for this trial
    trialIdx = trialOrder(iTrial);
    soa = p.soas(trials(trialIdx,1));
    targetPresent = p.targetPresAbs(trials(trialIdx,2));
    
    if soa < 0
        tex1 = surroundTex;
        if targetPresent == 1
            tex2 = targetTex;
        else
            tex2 = blankTex;
        end
        cueIm1SOA = p.cueTargetSOA + soa; % surround
        cueIm2SOA = p.cueTargetSOA; % target
    else
        if targetPresent == 1
            tex1 = targetTex;
        else
            tex1 = blankTex;
        end
        if soa==p.targOnlyCode % target only
            tex2 = blankTex;
            soa = 3*p.refRate; % set soa to a reasonable number
        else
            tex2 = surroundTex;
        end
        cueIm1SOA = p.cueTargetSOA; % target
        cueIm2SOA = p.cueTargetSOA + soa; % surround
    end
    
    % Select random jitter for stim 1 onset
    % ** consider using a flat hazard function instead **
    jitter = 0.2*rand; % uniform random delay between 0-200 ms
    
    % Present cue
    Screen('FillRect', win, [255 255 255], fixRect);
    timeCue = Screen('Flip', win);
    
    % Present images
    Screen('DrawTexture', win, tex1, [], imRect);
    Screen('FillRect', win, [255 255 255], fixRect);
    timeIm1 = Screen('Flip', win, timeCue + jitter + cueIm1SOA - slack);
    
    if abs(soa)>p.imageDur
        Screen('FillRect', win, white*p.backgroundColor);
        Screen('FillRect', win ,[255 255 255], fixRect);
        timeBlank1 = Screen('Flip', win, timeIm1 + p.imageDur - slack);
    end
    
    Screen('DrawTexture', win, tex2, [], imRect);
    Screen('FillRect', win, [255 255 255], fixRect);
    timeIm2 = Screen('Flip', win, timeCue + jitter + cueIm2SOA - slack);
    
    Screen('FillRect', win, white*p.backgroundColor);
    Screen('FillRect', win ,[255 255 255] ,fixRect);
    timeBlank2 = Screen('Flip', win, timeIm2 + p.imageDur - slack);
    
    Screen('FillRect', win, [0 0 255], fixRect);
    timeRespCue = Screen('Flip', win, timeBlank2 + p.postStimCushion - slack);
    
    % Find the actual times the target and surround were presented
    if soa < 0
        timeSurround = timeIm1;
        timeTarget = timeIm2;
    else
        timeTarget = timeIm1;
        timeSurround = timeIm2;
    end
    
    % Collect response
    [secs, keyCode] = KbWait(devNum);
    rt = secs - timeTarget;
    
    % Feedback
    responseKey = find(p.keyCodes==find(keyCode));
    response = p.targetPresAbs(responseKey);
    if response==targetPresent;
        correct = 1;
    else
        correct = 0;
    end
    
    % Store trial info
    trials(trialIdx,3) = soa;
    trials(trialIdx,4) = timeSurround - timeTarget;
    trials(trialIdx,5) = rt;
    trials(trialIdx,6) = responseKey;
    trials(trialIdx,7) = response;
    trials(trialIdx,8) = correct;
    trials(trialIdx,9) = jitter;
end
toc

% Show end of block feedback
% acc = mean(trials(:,6));

% Store expt info
expt.trialOrder = trialOrder;

results.trials_headers = trials_headers;
results.trials = trials;

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


