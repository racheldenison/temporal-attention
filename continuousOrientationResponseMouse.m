function [angle rt] = continuousOrientationResponseMouse(window, white, targetTex, targetRect, phRect, p)

% pick a starting orientation at random
angle = floor(180*rand);

% set mouse position and get buttons
startX = round(targetRect(1) + targetRect(3)/2);
startY = round(targetRect(2) + targetRect(4)/2);
SetMouse(startX, startY);
[x, y, buttons] = GetMouse(window);

% make sure left mouse button is not pressed to start
while buttons(1)~=0
    buttons(1) = 0;
end
WaitSecs(0.2); % clear the buffer
timeStart = GetSecs;
while buttons(1)==0
    xold=x;
    [x, y, buttons] = GetMouse;
    
    if xold-x>0 %&& buttons(1)==1
        angle=angle-1;
        if xold-x>5
            angle=angle-3;
            if xold-x>15
                angle=angle-6;
            end
        end
    elseif xold-x<0 %&& buttons(1)==1
        angle=angle+1;
        if xold-x<-5
            angle=angle+3;
            if xold-x<-15
                angle=angle+6;
            end
        end
    end
    
    DrawFormattedText(window, 'x', 'center', 'center', white);
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
    Screen('DrawTexture', window, targetTex, [], targetRect, angle);
    Screen('Flip', window);
end
timeEnd = GetSecs;
rt = timeEnd - timeStart;

angle = mod(angle,180);
WaitSecs(0.2);

