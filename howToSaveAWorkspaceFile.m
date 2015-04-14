% howToSaveAWorkspaceFile.m

% if you end a run early, save the data in TEMP (all the data collected to
% that point) as a sensibly named WORKSPACE file, subjectID_WORKSPACE.mat
%
% To run a session from a saved workspace, set the workspace variable in
% rd_temporalAttention.m

load data/TEMP.mat

workspaceFile = sprintf('data/%s_WORKSPACE.mat', subjectID);

save(workspaceFile)