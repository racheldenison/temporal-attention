function cost = rd_DoGCost(p, xdata, ydata)

% p is a vector: [mu sd amp b]

% mu = p(1);
% sd = p(2);
% amp = p(3);
% b = p(4);

sd = p(1);
amp = p(2);

mu = 0;
b = 0;

xgrid = -90:90;

for i=1:length(xdata)
    f = DoG(mu, sd, amp, b, xgrid);
    error(i) = ydata(i) - f(xgrid==xdata(i));
end

cost = sum(error.^2);