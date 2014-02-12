function imMasked = maskWithGaussian(im, sz, sigma)

g = make2DGaussianOval(sz, sz, round(sz/2), round(sz/2), ...
    sigma, sigma, 1);

imZeroCentered = im - mean(im(:));
imScaled = imZeroCentered.*g;
imMasked = imScaled + mean(im(:));