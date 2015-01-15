% texture_ssvef_test.m


%% setup
imSize = 1; % degrees
contrast = 1;
orientation = 0;
orientBandwidth = 10;
spatialFrequency = 4;
sfBandwidth = .5;
pixelsPerDegree = 100;
maskWithAperture = 1;

nIm = 100;

%% make images
for iIm = 1:nIm
    im1(:,:,iIm) = makeFilteredNoise(imSize, contrast, 0, orientBandwidth, ...
        spatialFrequency, sfBandwidth, pixelsPerDegree, maskWithAperture);
    im2(:,:,iIm) = makeFilteredNoise(imSize, contrast, 90, orientBandwidth, ...
        spatialFrequency, sfBandwidth, pixelsPerDegree, maskWithAperture);
end

c1 = (im1 + im2)./2;
c2 = ((1-im1) + (1-im2))./2;

%% make targets
t1 = makeFilteredNoise(imSize, contrast, -7, 3, ...
        spatialFrequency, sfBandwidth, pixelsPerDegree, maskWithAperture);
t2 = makeFilteredNoise(imSize, contrast, 97, 3, ...
        spatialFrequency, sfBandwidth, pixelsPerDegree, maskWithAperture);
    
%% play images
for iIm = 1:nIm
    if iIm == 30
        s = t1;
    elseif iIm == 70
        s = t2;
    else
        s = c1(:,:,iIm);
    end
    imshow(s)
    pause(0.05)
%     imshow(c2(:,:,iIm))
%     pause(0.05)
end


