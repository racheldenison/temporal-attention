function [targetRot rt startingAngle] = continuousOrientationResponseKeys(window, white, targetTex, targetRect, phRect, p)

% a couple more initializations
adjustments = [-1 1];

% pick a starting orientation at random
targetRots = randperm(180)-1;
targetRot = targetRots(1);
startingAngle = targetRot;

doneKeyPressed = 0;
timeStart = GetSecs;
while ~doneKeyPressed
    adjustment = 0;
    
    DrawFormattedText(window, 'x', 'center', 'center', white);
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    Screen('DrawTexture', window, targetTex, [], targetRect, targetRot);
    Screen('Flip', window);
    
    % check for orientation adjustment or final answer response
    [keyIsDown, seconds, keyCode] = KbCheck;
    if keyIsDown
        keyPressed = find(keyCode);
        doneKeyPressed = p.doneKey==keyPressed;
        adjustment = adjustments(p.keyCodes==keyPressed);
        if isempty(adjustment)
            adjustment = 0;
        end
    end
    
    % adjust orientation
    targetRot = targetRot + adjustment;
end
timeEnd = GetSecs;
rt = timeEnd - timeStart;

targetRot = mod(targetRot,180);

