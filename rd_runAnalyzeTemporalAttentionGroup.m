% rd_runAnalyzeTemporalAttentionGroup.m

% subjectInits = {'dg','sl','mr','ly','pv','ek','gk','md','ax'};
subjectInits = {'pv','ek','gk','md','ax'};
nSubjects = numel(subjectInits);

soas = {'100','300','800'};

for iSubject = 1:nSubjects
    subjectInit = subjectInits{iSubject};
    for iSOA = 1:numel(soas)
        soa = soas{iSOA};
        subject = sprintf('%s_cbD10_tilt*_tc16-64_soa1000-1%s', subjectInit, soa);
        rd_runAnalyzeTemporalAttention(subject);
    end
end