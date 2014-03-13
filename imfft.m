 function Y = imfft(X)

% Y = imfft(X)
%
% Performs 2D FFT on an image and rearranges result to place low frequencies centrally.

Y = fftshift(fft2(X));