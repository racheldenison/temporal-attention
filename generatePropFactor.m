function f = generatePropFactor(propA)
%
% f = generatePropFactor(propA)
%
% propA gives the proportion of A trials out of total trials
% f is the factor, which should have prop(A) 1s and prop(1-A) 2s

n = 100; % can approximate proportions to the nearest 1%

bigf = zeros(1,n);
bigf(1:round(propA*100)) = 1;

nA = nnz(bigf==1);
nB = nnz(bigf==0);

f = [ones(1,nA/gcd(nA,nB)) 2*ones(1,nB/gcd(nA,nB))];

if abs(nnz(f==1)/numel(f) - propA) > 0.001;
    warning('generatePropFactor can only approximate proportions to the nearest 1%. check propA.')
end
