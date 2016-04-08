function bullseye = buildBullseye(spatialFrequency,sizeDegrees,pixelsPerDegree,maskContrast)

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

[x y] = meshgrid(dimArray{2}, dimArray{1});

circles =  sin(radiansPerPixel*sqrt(y.^2 + x.^2) + pi/2);

bullseye = (circles*maskContrast)/2 + .5;

