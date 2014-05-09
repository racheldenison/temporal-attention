% rd_blackWhite.m
%
% This is designed to be a very simple "experiment" that can be used to
% test and/or demonstrate the use of rd_eyeLink.m.

subjectID = 'ellipse';
eyeDataDir = 'eyedata';
eyeTracking = 1;

% eyeFile = sprintf('%s%s', subjectID, datestr(now, 'mmdd'));
eyeFile = subjectID;

nTrials = 10;
rad = 70; % radius of allowable eye movement in pixels
stimDur = 2; % duration of black and white screens

%% Screen
screenNumber = max(Screen('Screens'));
[window rect] = Screen('OpenWindow', screenNumber);
% [window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]); % for testing
[cx cy] = RectCenter(rect);
Screen('TextSize', window, 24);
Screen('TextColor', window, 255);
Screen('TextFont', window, 'Verdana');
Screen('FillRect', window, 128)
Screen('Flip', window);

%% Initialize eye tracker
if eyeTracking
    [el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
    if exitFlag
        return
    end
    
    %% Calibrate eye tracker
    [cal exitFlag] = rd_eyeLink('calibrate', window, el);
    if exitFlag
        return
    end
end

%% Present gray screen
Screen('FillRect', window, 128)
DrawFormattedText(window, '+', 'center', 'center', [255 0 0]);
Screen('Flip', window);

%% Present trials
for iTrial = 1:nTrials
    Beeper('low')
    KbWait;
    
    fprintf('\n\nTrial %d', iTrial)
    % present fixation
    DrawFormattedText(window, '+', 'center', 'center', [255 0 0]);
    timeFix = Screen('Flip', window);
    
    % start the trial when the eyetracker is recording and the subject is
    % holding fixation
    if eyeTracking
        rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});
    end
    
    % present first stimulus: BLACK
    Screen('FillRect', window, 0)
    DrawFormattedText(window, '+', 'center', 'center', [255 0 0]);
    %     DrawFormattedText(window, 'STIM 1', cx-200, 'center');
    timeStim1(iTrial) = Screen('Flip', window);
    if eyeTracking
        Eyelink('Message', 'EVENT_BLACK');
    end
    
    stopThisTrial = 0;
    while GetSecs < timeStim1(iTrial) + stimDur && ~stopThisTrial
        WaitSecs(.01);
        if eyeTracking
            fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            if ~fixation
                fprintf('\nBROKE FIXATION! (stim 1)')
                Beeper('high')
                stopThisTrial = 1;
            end
        end
    end
    
    if stopThisTrial
        continue
    end
    
    % present second stimulus: WHITE
    Screen('FillRect', window, 255)
    DrawFormattedText(window, '+', 'center', 'center', [255 0 0]);
    %     DrawFormattedText(window, 'STIM 2', cx+200, 'center');
    timeStim2(iTrial) = Screen('Flip', window);
    if eyeTracking
        Eyelink('Message', 'EVENT_WHITE');
    end
    
    while GetSecs < timeStim2(iTrial) + stimDur && ~stopThisTrial
        WaitSecs(.01);
        if eyeTracking
            fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            if ~fixation
                fprintf('\nBROKE FIXATION! (stim 2)')
                Beeper('high')
                stopThisTrial = 1;
            end
        end
    end
    
    if eyeTracking
        Eyelink('Message', 'TRIAL_END');
    end
end

%% Final blank
Screen('FillRect', window, 128)
Screen('Flip', window);

%% Save the eye data and shut down the eye tracker
if eyeTracking
    if ~exist(eyeDataDir,'dir')
        mkdir(eyeDataDir)
    end
    rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});
end

%% Close screen
Screen('CloseAll')

