% gaussian_notes.m

sz = 100;
sigma = 20;

g = make2DGaussianOval(sz, sz, round(sz/2), round(sz/2), ...
    sigma, sigma, 1);
g50 = g(:,50);

x = 1:100;
y = normpdf(x,50,20);
y1 = normpdf(x,50,20*sqrt(2));

figure 
hold on
plot(g50)
plot(y,'r')
plot(y*(1/max(y)),'r')
plot(y1*(1/max(y1)),'g')