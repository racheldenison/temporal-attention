% rd_writeTemporalAttentionDataFile.m

%% multi SOA
% D = load('data/E2_SOA_cbD6_run98_N5_workspace_20180731.mat');
D = load('data/E2_SOA_cbD6_run98_N5_workspace_20191003.mat');
data = D.dpData; % dpData, rtData
subjects = D.subjectInits;
soas = D.t1t2soa;
ts = D.intervalNames;
vs = D.cueNames;

fileID = fopen('data/E2_SOA_cbD6_run98_N5_dprime.txt','w');
fprintf(fileID,'%s %s %s %s %s\n','subject','target','soa','validity','dprime');

for iS = 1:numel(subjects)
    subject = subjects{iS}; % initials
%     subject = sprintf('s%d', iS); % anonymized subject number
    for iT = 1:numel(ts)
        t = ts{iT};
        for iSOA = 1:numel(soas)
            soa = soas(iSOA);
            for iV = 1:numel(vs)
                v = vs{iV};
                val = data{iT}(iV,iSOA,iS);
                fprintf(fileID,'%s %s %d %s %1.4f\n', ...
                    subject, t, soa, v, val);
            end
        end
    end
end
    
fclose(fileID);

%% car format
strs = {};
vals = [];
idx = 1;
for iT = 1:numel(ts)
    t = ts{iT};
    for iSOA = 1:numel(soas)
        soa = soas(iSOA);
        for iV = 1:numel(vs)
            v = vs{iV};
            strs{idx} = sprintf('%ssoa%d%s',t,soa,v);
            vals(:,idx) = data{iT}(iV,iSOA,:);
            idx = idx + 1;
        end
    end
end

% data matrix
fileID = fopen('data/E2_SOA_cbD6_run98_N5_rt_cardata.txt','w');
fprintf(fileID,'subject ');
for iStr = 1:numel(strs)
    fprintf(fileID,'%s ',strs{iStr});
end

for iS = 1:numel(subjects)
    fprintf(fileID,'\n%s ',subjects{iS});
    for iStr = 1:numel(strs)
        fprintf(fileID,'%1.4f ',vals(iS,iStr));
    end
end

% idata matrix
fileID = fopen('data/E2_SOA_cbD6_run98_N5_rt_caridata.txt','w');
fprintf(fileID, '%s %s %s\n', 'target', 'soa', 'validity'); 
for iStr = 1:numel(strs)
    str = strs{iStr};
    fprintf(fileID, '%s %s %s\n', str(1:2), str(3:8), str(9:end));
end
