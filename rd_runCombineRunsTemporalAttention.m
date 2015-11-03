% rd_runCombineRunsTemporalAttention.m

% subjectInits = {'dg','sl','mr','ly','pv','ek','gk','md','ax'};
% soaStrs = {'1100','1300','1800'};
subjectInits = {'rd'};
soaStrs = {'1450'};

for iSubject = 1:numel(subjectInits)
    subjectInit = subjectInits{iSubject};
    
    for iSOA = 1:numel(soaStrs)
        soaStr = soaStrs{iSOA};
    
        subject = sprintf('%s_cbD6*_tc64_soa1000-%s', subjectInit, soaStr);

        rd_combineRunsTemporalAttention(subject);
        close all
    end
end