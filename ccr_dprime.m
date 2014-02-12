function [dprime,criterion] = ccr_dprime(h,fa,option)

if (h >= 1) | (fa >= 1) | (h <= 0) | (fa <= 0)
    error('h and fa must be non-ceiling proportions, i.e. >= 0 & =< 1!')
end

if isempty(option) | ~ischar(option)
    option = input('what do you want? ''yesno'' or ''2afc'': ');
end

zh = norminv(h,0,1); zfa = norminv(fa,0,1);

switch option
    case 'yesno'
        dprime = zh - zfa;
        criterion = -0.5*(zh+zfa);
    case '2afc'
        dprime = (zh - zfa)/sqrt(2);  
        criterion = -0.5*(zh+zfa);
    otherwise
        error('option not recognised! please use ''yesno'' or ''2afc''!');
end

return