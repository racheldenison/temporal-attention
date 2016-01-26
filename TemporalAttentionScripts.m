% TemporalAttentionScripts

rd_checkParamsTemporalAttention
rd_combineRunsTemporalAttention

% to re-run analysis
% eg. look at same/diff axis, clean RT
rd_runAnalyzeTemporalAttentionGroup
rd_runAnalyzeTemporalAttention

% group analysis for E0
rd_analyzeTemporalAttentionGroup

% multi-SOA
rd_plotTemporalAttentionMultiSOAGroup
% -> settings in:
rd_plotTemporalAttentionMultiSOA 
    % can reanalyze with different T1T2Axis here

% also, for contrast on x-axis
rd_plotTemporalAttentionMultiSOAGroupCRF
rd_plotTemporalAttentionMultiSOACRF 

% multi-SOA resampling
rd_plotTemporalAttentionMultiSOABootstrapErrorBars
rd_plotTemporalAttentionMultiSOAResampledNull
rd_resampleTemporalAttentionMultiSOA

% main analysis
rd_analyzeTemporalAttention
    % can perform an extra selection step here (e.g. same/diff orientation)
    


    