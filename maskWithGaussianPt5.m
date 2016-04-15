function imMasked = maskWithGaussianPt5(im, sz, sigma)

g = make2DGaussianOval(sz, sz, round(sz/2), round(sz/2), ...
    sigma, sigma, 1);

imZeroCentered = im - .5;
imScaled = imZeroCentered.*g;
imMasked = imScaled + .5;