function Y = imifft(X)

% Y = imifft(X)
%
%imifft.m - Inverse imFFT.

Y = abs(ifft2(ifftshift(X)));