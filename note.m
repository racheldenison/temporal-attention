function y = note(f1,fs,seconds)

% Here's an alternate version which allows multiple frequency components
% for chords, harmonics, etc. For example to create an A major chord you 
% pass the vector [220 277.12 329.63 440] as the f1 argument.
 
%This determines how many frequency components the sound will have
len = length(f1);
 
%This piece of code sets up a vector that, when cycled through at the sampling frequency (fs), will 
%take "seconds" seconds from start to finish
n = 1:seconds*fs;
 
%This creates the vector y
y(n) = 0;
for i = 1:len
  y = y + sin(2*pi*(f1(i)/fs)*n);
end

