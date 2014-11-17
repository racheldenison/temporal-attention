% fixRespErrorGr90.m

fileName = '/E3_adjust/bl/bl_a1_tc100_soa1000-1250_run02_TemporalAttentionAdjust_T1T2all_20141114';
load([pathToExpt('data') fileName])

responseError = expt.trials(:,11);
responseError1 = responseError;
responseError1(responseError<-90) = ...
    responseError1(responseError<-90)+180;
respErrorAbs1 = abs(responseError1);
expt.trials(:,11:12) = [responseError1 respErrorAbs1];

save([pathToExpt('data') fileName '_fixed'], 'expt', 'results')

