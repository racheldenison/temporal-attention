% howToTransferEyeData.m

% Ctrl-C to quit experiment program
% DO NOT clear all
% run this program

load('data/TEMP.mat')

% eyeDataDir = 'eyedata';
% eyeFile = 'testeyefile.edf';

Eyelink('StopRecording');
Eyelink('ReceiveFile', eyeFile, eyeDataDir, 1);
Eyelink('CloseFile');
Eyelink('Shutdown');

% rename eye file
eyeFileFull = sprintf('%s/%s_PART1_TemporalAttention_%s.edf', eyeDataDir, subjectID, datestr(now, 'yyyymmdd'));
copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)