% rd_eyeTrackTester.m
%
% This is designed to be a very simple "experiment" that can be used to
% test and/or demonstrate the use of rd_eyeLink.m.

subjectID = 'test';
eyeDataDir = 'eyedata';

eyeFile = sprintf('%s%s', subjectID, datestr(now, 'mmdd'));

nTrials = 10;
rad = 40; %%% pixels? % radius of allowable eye movement

%% Screen
screenNumber = max(Screen('Screens'));
[window rect] = Screen('OpenWindow', screenNumber);
% [window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]); % for testing
[cx cy] = RectCenter(rect);
Screen('TextSize', window, 24);
Screen('TextColor', window, 255);
Screen('TextFont', window, 'Verdana');

%% Initialize eye tracker
[el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
if exitFlag
    return
end

%% Calibrate eye tracker
[cal exitFlag] = rd_eyeLink('calibrate', window, el);
if exitFlag
    return
end

%% Present trials
for iTrial = 1:nTrials
    % wait for a keypress, then go on to the next trial
    KbWait;
    
    fprintf('\n\nTrial %d', iTrial)
    % present fixation
    DrawFormattedText(window, '+', 'center', 'center');
    timeFix = Screen('Flip', window);
    
    % start eye recording for this trial
    % don't start the trial until the subject is holding fixation
    rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});
    
    % present first stimulus
    DrawFormattedText(window, '+', 'center', 'center');
    DrawFormattedText(window, 'STIM 1', cx-200, 'center');
    timeStim1(iTrial) = Screen('Flip', window);
    
    stopThisTrial = 0;
    while GetSecs < timeStim1(iTrial) + 0.45 && ~stopThisTrial
        WaitSecs(.01);
        fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
        if ~fixation
            % do fixation break tasks:
            % this could include adding the broken trial to the end
            fprintf('\nBROKE FIXATION!\n')
            stopThisTrial = 1;
        end
    end
    fix1(iTrial) = fixation;
    % alt: fixationBreakTasks(fixation)
    
    if stopThisTrial
        continue
    end
    
    % present second stimulus
    DrawFormattedText(window, '+', 'center', 'center');
    DrawFormattedText(window, 'STIM 2', cx+200, 'center');
    timeStim2(iTrial) = Screen('Flip', window, timeStim1(iTrial) + 0.5);
    
    while GetSecs < timeStim2(iTrial) + 0.45 && ~stopThisTrial
        WaitSecs(.01);
        fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
        if ~fixation
            fprintf('\n\nBROKE FIXATION!\n')
            stopThisTrial = 1;
        end
    end
    fix2(iTrial) = fixation;
    
    % stop eye recording for this trial
    rd_eyeLink('trialstop', window);
end

 %% Save the eye data and shut down the eye tracker
if ~exist(eyeDataDir,'dir')
    mkdir(eyeDataDir)
end
rd_eyeLink(window, 'eyestop', {eyeFile, eyeDataDir});

%% Close screen
Screen('CloseAll')

