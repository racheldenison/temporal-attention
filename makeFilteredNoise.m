function filteredNoiseIm = makeFilteredNoise(imSize, contrast, ...
    orientation, orientBandwidth, ...
    spatialFrequency, sfBandwidth, ...
    pixelsPerDegree, maskWithAperture)

% function filteredNoiseIm = makeFilteredNoise(imSize, contrast, ...
%     orientation, orientBandwidth, ...
%     spatialFrequency, sfBandwidth, ...
%     pixelsPerDegree, maskWithAperture)
%
% INPUTS:
% imSize is the 1-element image size in degrees
% orientation is the main orientation in degrees
% orientBandwidth is the width of orientation bandwidth in degrees
% contrast is the contrast from 0-1
% spatialFrequency is the main spatial frequency in cpd
% sfBandwidth is the width of the spatial frequency bandwidth in cpd
% pixelsPerDegree is pixels per degree of visual angle
% maskWithAperture is 1 if want to mask with an aperture, 0 if not. Note
%   that the imSize will always determine the size of the visible noise
%   patch. If we mask with an aperture, the size of the entire image will be
%   1.3 times larger than if we don't (the aperture is added around the noise patch).

%% example inputs
if nargin==0
    imSize = 1; % degrees
    contrast = 1;
    orientation = 0;
    orientBandwidth = 10;
    spatialFrequency = 4;
    sfBandwidth = .5;
    pixelsPerDegree = 100;
    maskWithAperture = 1;
end

% viewDist = 60;
% screenPx = 1024;
% screenDim = 20;

%% minor calculations
% pixelsPerDegree = ang2pix(1, screenDim, screenPx, viewDist, 'central');
% dimPerDegree = pixelsPerDegree*(screenDim/screenPx);
% samplingRate = pixelsPerDegree/dimPerDegree;
samplingRate = pixelsPerDegree;
fNyquist = samplingRate/2;
fLow = (spatialFrequency-sfBandwidth/2)/fNyquist;
fHigh = (spatialFrequency+sfBandwidth/2)/fNyquist;
sz = round(imSize*pixelsPerDegree);

%% make smoothing filter
[x1,y1]=meshgrid(-10:10,-10:10);
sigma=1/3*6;
smoothfilter = 1*exp((-2.77*x1.^2)/(2.35*sigma)^2).*exp((-2.77*y1.^2)/(2.35*sigma)^2); 

%% make aperture
if maskWithAperture
    ap = ones(round(sz * 1.3), round(sz * 1.3));
    center = sz * 1.3/2;
    for i=1:sz * 1.3
        for j=1:sz * 1.3
            R(i,j)=sqrt(((i-center).^2)+((j-center).^2));
        end
    end
    out_index_L=find(R>sz/2);
    ap(out_index_L) = 0;
    aperture = filter2(fspecial('gaussian', round(sz/2), round(sz/8)), ap);
else
    aperture = ones(sz);
end

%% make filters
oFilter = OrientationBandpass(length(aperture), orientation - orientBandwidth/2, orientation + orientBandwidth/2);
fFilter = Bandpass2(length(aperture), fLow, fHigh);

oFilter = filter2(smoothfilter, oFilter);
fFilter = filter2(smoothfilter, fFilter);

noise = normrnd(0, 1, length(aperture), length(aperture));
fn = fftshift(fft2(noise));
filteredNoise = real(ifft2(ifftshift(oFilter.*fFilter.*fn)));

filteredNoiseS = Scale(filteredNoise);

filteredNoiseIm = contrast * (filteredNoiseS-0.5) .* aperture + 0.5;

%% show the image
% figure
% subplot(2,1,1)
% imshow(filteredNoiseS)
% subplot(2,1,2)
% imshow(filteredNoiseIm)

        