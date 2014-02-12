function grating = buildColorGrating(pixelsPerDegree, sizeDegrees, ...
    spatialFrequency, tiltInDegrees, phase, michelsonContrast, ...
    squareWave, colorOption, redWeight, greenWeight)

% function grating = buildColorGrating(pixelsPerDegree, ...
%     spatialFrequency, tiltInDegrees, phase, michelsonContrast, ...
%     squareWave, colorOption, redWeight, greenWeight)
%
% Builds a black/white or red/green grating
%
% Inputs:
% pixelsPerDegree [=99]: pixels per degree of visual angle
% sizeDegrees [=[2 2]]: size [height width] in degrees of visual angle
% spatialFrequency [=3]
% tiltInDegrees [=0]
% phase [=0]
% michelsonContrast [=1]: contrast of the grating, 0-1
% squareWave [=0]: 1 for square wave gratings, 0 for sine wave gratings
% colorOption [='bw']: 'bw' for black/white or 'rg' for red/green
% redWeight [=1]: 0-1
% greenWeight [=1]: 0-1
%
% Outputs:
% grating: 2D black/white or 3D red/green grating with range 0-1 for each
% color channel
%
% Calling the function with zero arguments will generate 100% contrast
% black/white grating with default parameter values.
%
% squareWave, colorOption, redWeight, and greenWeight are all optional, 
% but if you specify one of them, you must also specify all the previous 
% ones.

if nargin==0
    pixelsPerDegree = 36.5; % 99; % this is with NEC Monitor with subject sitting 5 feet from the screen. 1280 x 1024 pixel fullscreen.
    sizeDegrees = [2 3];
    spatialFrequency = .5;
    tiltInDegrees = 0; 
    phase = 0;
    michelsonContrast = 1;
end
if nargin < 7
    squareWave = 0;
end
if nargin < 8
    colorOption = 'rg'; 
end
if nargin < 9
    redWeight = 1;
    greenWeight = 1;
end

pixelsPerCycle = pixelsPerDegree / spatialFrequency; % How many pixels will each period/cycle occupy?
cyclesPerPixel = 1/pixelsPerCycle; % How many periods/cycles are there in a pixel?
radiansPerPixel = cyclesPerPixel * (2 * pi); % = (periods per pixel) * (2 pi radians per period)

% Set diameter of grating
sizePixels = sizeDegrees * pixelsPerDegree;

% Size of grid
for iDim = 1:numel(sizePixels)
    dimOfGrid = sizePixels(iDim);   %the next lines make sure that it is a whole, even number so matrix indices don't choke.
    dimOfGrid = round (dimOfGrid);
    if mod (dimOfGrid, 2) ~= 0
        dimOfGrid = dimOfGrid + 1 ;
    end
    
    halfDimOfGrid =  (dimOfGrid / 2);
    dimArray{iDim} = (-halfDimOfGrid) : halfDimOfGrid;  % widthArray is used in creating the meshgrid.
    sizeOfGrid(iDim) = length(dimArray{iDim});
end

% Set tilt (angle of grating)
tiltInRadians = tiltInDegrees * pi / 180; % The tilt of the grating in radians.

% ---------- Image Setup ----------
% Creates a two-dimensional square grid.  For each element i = i(x0, y0) of
% the grid, x = x(x0, y0) corresponds to the x-coordinate of element "i"
% and y = y(x0, y0) corresponds to the y-coordinate of element "i"
[x y] = meshgrid(dimArray{2}, dimArray{1});

% if we want dot in middle: ...|((x.^2 + y.^2) <  (pixelsPerDegree *.05)^2 ) ;

% ---------- Build grating -----------
a1 = cos(tiltInRadians) * radiansPerPixel;
b1 = sin(tiltInRadians) * radiansPerPixel;

sinwav = sin(a1*x+b1*y+phase);
if squareWave
    sinwav(sinwav<=0) = -1;
    sinwav(sinwav>0) = 1;
end

imageMatrix = .5 + .5*(michelsonContrast * sinwav);

switch colorOption
    case 'bw'
        grating = imageMatrix;  
    case 'rg'
        imageMatrixRGB(:,:,1) = (imageMatrix * redWeight);   %set red gun on
        imageMatrixRGB(:,:,2) = (1 - imageMatrix) * greenWeight; % set green to negative spaces
        imageMatrixRGB(:,:,3) = zeros(sizeOfGrid); % ignore blue
        
        grating = imageMatrixRGB;
end

% imshow(grating)

