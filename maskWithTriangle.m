function imMasked = maskWithTriangle(im, x0, y0, triW, triH, blurSize)

t0 = makeUprightTriangle(size(im,2), size(im,1), ...
    x0, y0, triW, triH);

t = blurMaskEdges(t0, blurSize);

imZeroCentered = im - mean(im(:));
imScaled = imZeroCentered.*t;
imMasked = imScaled + mean(im(:));