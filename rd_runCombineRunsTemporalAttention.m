% rd_runCombineRunsTemporalAttention.m

subjectInits = {'dg','sl','mr','ly','pv','ek','gk','md','ax'};
soaStrs = {'1100','1300','1800'};
% subjectInits = {'ax'};

for iSubject = 1:numel(subjectInits)
    subjectInit = subjectInits{iSubject};
    
    for iSOA = 1:numel(soaStrs)
        soaStr = soaStrs{iSOA};
    
        subject = sprintf('%s_cbD10_tilt*_tc16-64_soa1000-%s', subjectInit, soaStr);

        rd_combineRunsTemporalAttention(subject);
        close all
    end
end