function [out, exitFlag] = rd_eyeLink(window, command, in)

% Possible commands and their ins & outs:
% 'eyestart'
%   in = eyeFile
%   out = el
%
% 'calibrate'
%   in = el
%   out = cal
%
% 'trialstart'
%   in = {el, trialNum, cx, cy, rad}
%   out = []
%
% 'fixholdcheck'
%   in = {cx, cy, rad}
%   out = fixation
%
% 'fixcheck'
%   in = {cx, cy, rad}
%   out = fixation
%
% 'driftcorrection'
%   in = {el, cx, cy}
%   out = driftCorrection
%
% 'trialstop'
%   in = []
%   out = []
%
% 'eyestop'
%   in = {eyeFile, eyeDataDir}
%   out = []

%% Initializations
% assume no output unless some is given
out = [];

% assume everything goes ok (exitFlag=0) until proven otherwise        
exitFlag = 0;
        
%% Do command
switch command
    case 'eyestart'
        %% start eyetracker
        eyeFile = in;
        
        % First check if we can get a connection
        if EyelinkInit()~= 1
            fprintf('Couldn''t initialize connection with eyetracker! Exiting ...\n');
            return
        end
        
        % Set up the eyetracker
        el = EyelinkInitDefaults(window);
        el.backgroundcolour = [128 128 128]; %%% needed?
        
        Eyelink('Command', 'file_sample_data = LEFT,RIGHT,GAZE,AREA');
        % Eyelink('Command', 'calibration_type = HV9');
        [v vs] = Eyelink('GetTrackerVersion');
        fprintf('\nRunning experiment on a ''%s'' tracker.', vs );
        
        % Start the eye file
        edfFile = sprintf('%s.edf', eyeFile);
        edfFileStatus = Eyelink('Openfile', edfFile);
        if edfFileStatus==0
            fprintf('\nEye file opened ok.\n\n')
        else
            fprintf('Cannot open .edf file. Exiting ...');
            Screen('CloseAll')
            exitFlag = 1;
            return
        end
        
        out = el; % return the el structure as output
        
    case 'calibrate'
        %% calibrate eyetracker
        el = in;
        
        calString = sprintf('Eye tracker calibration:\n\nPlease fixate the center of the dot!\n\nPress ''space'' to start or ''q'' to quit!');
        DrawFormattedText(window, calString, 'center', 'center', 1, []);
        Screen('Flip', window, 0, 1);
        
        contKey = '';
        while isempty(find(strcmp(contKey,'space'), 1))
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, keyTime, keyCode] = KbCheck;
            end
            contKey = KbName(find(keyCode));
        end
        if strcmp(contKey,'q')
            ListenChar(0);
            ShowCursor;
            Screen('CloseAll')
            fclose('all');
            fprintf('User ended program');
            exitFlag = 1;
            return
        end
        Screen('Flip', window, 0, 1);
        
        cal = EyelinkDoTrackerSetup(el);
        if cal==el.TERMINATE_KEY
            exitFlag = 1;
            return
        end
        
        out = cal;
        
    case 'trialstart'
        %% trial start
        % start recording and wait until the subject is fixating
        el = in{1};
        trialNum = in{2};
        cx = in{3};
        cy = in{4};
        rad = in{5};
        
        % Displays a title at the bottom of the eye tracker display
        Eyelink('command', 'record_status_message ''Starting trial %d''', trialNum);
        
        fixation = 0; % fixation ok?
        record = 0; % are we recording?
        recalibrate = 0; % do we need to recalibrate?
        
        % Start the trial only when 1) eyetracker is recording, 2) subject
        % is fixating
        while ~record && ~fixation
            % Start recording and check that it really worked
            while ~record
                Eyelink('StartRecording');	% start recording
                % start recording 100 msec before just to be safe
                WaitSecs(0.1);
                key=1;
                while key~=0
                    key = EyelinkGetKey(el); % dump any pending local keys
                end
                
                err=Eyelink('CheckRecording'); 	% check recording status
                if err==0
                    record = 1;
                    Eyelink('message', 'RECORD_START');
                else
                    record = 0;	% results in repetition of fixation check
                    Eyelink('message', 'RECORD_FAILURE');
                end
            end
            
            % Verify that the subject is holding fixation for some set
            % time before allowing the trial to start. A
            % timeout period is built into this function.
            % (woah, recursion)
            fixation = rd_eyeLink(window, 'fixholdcheck', {cx, cy, rad});
            
            %%% what are the conditions for recalibration?? %%%
            % Recalibrate if needed
            if recalibrate
                cal = EyelinkDoTrackerSetup(el);
                if cal==el.TERMINATE_KEY
                    exitFlag = 1;
                    return
                end
                record = 0; % start over
            end
        end
        
        Eyelink('Message', 'TRIAL_START %d', trialNum);
        Eyelink('Message', 'SYNCTIME');		% zero-plot time for EDFVIEW
        
    case 'fixholdcheck'
        %% check that fixation is held for some amount of time
        cx = in{1}; % x coordinate of screen center
        cy = in{2}; % y coordinate of screen center
        rad = in{3}; % acceptable fixation radius %%% in px?
        
        timeout = 1.00; % maximum fixation check time
        tFixMin = 0.20; % minimum correct fixation time
        
        % determine recorded eye
        evt = Eyelink('newestfloatsample');
        domEye = find(evt.gx ~= -32768);
        
        Eyelink('Message', 'EVENT_FixationCheck');
        
        tstart = GetSecs;
        fixation = 0; % is the subject fixating now?
        fixStart = 0; % has a fixation already started?
        tFix = 0; % how long has the current fixation lasted so far?
        
        t = tstart;
        while ((t-tstart) < timeout && tFix<=tFixMin)
            % get eye position
            evt = Eyelink('newestfloatsample');
            x = evt.gx(domEye);
            y = evt.gy(domEye);
            
            % check eye position
            if sqrt((x-cx)^2+(y-cy)^2)<rad
                fixation = 1;
            else
                fixation = 0;
            end
            
            % update duration of current fixation
            if fixation==1 && fixStart==0
                tFixStart = GetSecs;
                fixStart = 1;
            elseif fixation==1 && fixStart==1
                tFix = GetSecs-tFixStart;
            else
                fixStart = 0;
            end
            
            t = GetSecs;
        end
        
        out = fixation;
        
    case 'fixcheck'
        %% check fixation at one point in time
        cx = in{1}; % x coordinate of screen center
        cy = in{2}; % y coordinate of screen center
        rad = in{3}; % acceptable fixation radius %%% in px?
        
        % determine recorded eye
        evt = Eyelink('newestfloatsample');
        domEye = find(evt.gx ~= -32768);
        
        Eyelink('Message', 'EVENT_FixationCheck');
        
        % get eye position
        x = evt.gx(domEye);
        y = evt.gy(domEye);
        
        % check eye position
        if sqrt((x-cx)^2+(y-cy)^2)<rad
            fixation = 1;
        else
            fixation = 0;
            Eyelink('Message', sprintf('BROKE FIXATION'));
        end
        
        out = fixation;
        
    case 'driftcorrection'
        %% do a drift correction
        el = in{1};
        cx = in{2};
        cy = in{3};
        
        driftCorrection = EyelinkDoDriftCorrect(el, cx, cy, 1, 1);
        
        out = driftCorrection;
        
    case 'trialstop'
        %% stop recording
        Eyelink('StopRecording');
        Eyelink('Message','EVENT_trialEnd');
        
    case 'eyestop'
        %% get the eye file and close down the eye tracker
        eyeFile = in{1};
        eyeDataDir = in{2};
        
        Eyelink('ReceiveFile', eyeFile, eyeDataDir, 1); 
        Eyelink('CloseFile'); 
        Eyelink('Shutdown');
        
    otherwise
        error('[rd_eyeLink]: ''command'' argument not recognized. See help for available commands.')
end
