function y = rd_wrappednormpdf(x, m, sigma, n)

if nargin==0
    x = linspace(-pi,pi,100);
    m = 0;
    sigma = pi/4;
    n = 2;
end

a = zeros(size(x));
for i = -1000:1000
    a = a + exp((-abs(x - m - 2*pi*i).^n)/(2*sigma.^n));
end

y = 1/(sigma*sqrt(2*pi)).*a; 

[xsorted idx] = sort(x);
ysorted = y(idx);

areay = sum(ysorted(1:end-1).*abs(diff(xsorted)));
y = y./areay;
