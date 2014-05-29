function [dprime,criterion] = rd_dprime(h,fa,option,ceilopt)

if nargin < 4 || isempty(ceilopt)
    ceilopt = 'noadjust';
end

switch ceilopt
    case 'noadjust'
        if any(h >= 1) | any(fa >= 1) | any(h <= 0) | any(fa <= 0)
            error('h and fa must be non-ceiling proportions, i.e. >= 0 & =< 1!')
        end
    case 'adjust'
        h(h==1) = .99;
        h(h==0) = .01;
        fa(fa==1) = .99;
        fa(fa==0) = .01;
    otherwise
        error('ceilopt not recognized.')
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