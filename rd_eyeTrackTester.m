% rd_eyeTrackTester.m
%
% This is designed to be a very simple "experiment" that can be used to
% test and/or illustrate the use of rd_eyeLink.m.

subjectID = 'eyetest';
eyeDataDir = 'eyedata';

% eyeFile = sprintf('%s_%s', subjectID, datestr(now, 'yyyymmdd'));
eyeFile = 'eyetest';

nTrials = 10;
rad = 50; %%% pixels? % radius of allowable eye movement

%% Screen
screenNumber = max(Screen('Screens'));
[window rect] = Screen('OpenWindow', screenNumber);
% [window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]); % for testing
[cx cy] = RectCenter(rect);
Screen('TextSize', window, 24);
Screen('TextColor', window, 255);
Screen('TextFont', window, 'Verdana');

%% Initialize eye tracker
[el exitFlag] = rd_eyeLink(window, 'eyestart', eyeFile);
if exitFlag
    return
end

%% Calibrate eye tracker
[cal exitFlag] = rd_eyeLink(window, 'calibrate', el);
if exitFlag
    return
end

%% Present trials
for iTrial = 1:nTrials
    % present fixation
    DrawFormattedText(window, '+', 'center', 'center');
    timeFix = Screen('Flip', window);
    
    % start eye recording for this trial
    % don't start the trial until the subject is holding fixation
    rd_eyeLink(window, 'trialstart', {el, iTrial, cx, cy, rad});
    
    % present first stimulus
    DrawFormattedText(window, 'STIM 1', cx-200, 'center');
    timeStim1 = Screen('Flip', window);
    
    % check fixation right after stimulus presentation
    fixation = rd_eyeLink(window, 'fixcheck', {cx, cy, rad});
    if ~fixation
        % do fixation break tasks:
        % this could include adding the broken trial to the end
        % of the trial list, but here we will just add an extra trial
        nTrials = nTrials+1;
        continue % stop this trial and go on to the next one
    end
    % alt: fixationBreakTasks(fixation)
    
    % present second stimulus
    DrawFormattedText(window, 'STIM 2', cx+200, 'center');
    timeStim2 = Screen('Flip', window, timeStim1 + 0.25);
    
    % check fixation right after stimulus presentation
    fixation = rd_eyeLink(window, 'fixcheck', {cx, cy, rad});
    if ~fixation
        nTrials = nTrials+1;
        continue
    end
    
    % stop eye recording for this trial
    rd_eyeLink(window, 'trialstop');
    
    % wait for a keypress, then go on to the next trial
    KbWait;
end

%% Save the eye data and shut down the eye tracker
if ~exist(eyeDataDir,'dir')
    mkdir(eyeDataDir)
end
rd_eyeLink(window, 'eyestop', {eyeFile, eyeDataDir});


