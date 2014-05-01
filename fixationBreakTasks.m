function [stopThisTrial trialOrder, nTrials] = fixationBreakTasks(...
    fixation, window, fillColor, trialOrder, iTrial, nTrials)

fixBreakSound = soundFreqSweep(100,500,.1221);
if fixation==0
    stopThisTrial = 1;
    
    % blank the screen, make a sound, and give a time out
    soundsc(fixBreakSound)
    Screen('FillRect', window, fillColor);
    Screen('Flip', window);
    WaitSecs(1);
    
    % redo this trial at the end of the experiment
    % this can be easily done by appending the trial number to the end of
    % trialOrder
    trialOrder(end+1) = trialOrder(iTrial);
    nTrials = nTrials + 1;
else
    stopThisTrial = 0;
end
