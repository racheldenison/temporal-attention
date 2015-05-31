function cost = rd_DoGCost(p, xdata, ydata)

% p is a vector: [amp sd mu b]

switch numel(p)
    case 1
        amp = p(1);
        sd = 44;
        mu = 0;
        b = 0;
    case 2
        amp = p(1);
        sd = p(2);
        mu = 0;
        b = 0;
    case 4
        amp = p(1);
        sd = p(2);
        mu = p(3);
        b = p(4);
    otherwise
        error('p has wrong number of elements')
end

xgrid = -90:90;

for i=1:length(xdata)
    f = DoG(mu, sd, amp, b, xgrid);
    err(i) = ydata(i) - f(xgrid==xdata(i));
end

cost = sum(err.^2);