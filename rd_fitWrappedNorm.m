% rd_fitWrappedNorm.m

%% standard
% function handle
fn = @(data, mu, sigma, n) rd_wrappednormpdf(data, mu, sigma, n);
% fn = @(data, mu, sigma, n) wrappednormnpdf(data, mu, sigma, n);

% data to fit
x = normrnd(pi/8,pi/4,1,1000);
% x = normrnd(5,30,1,1000);

% mle
[phat, pci] = mle(x,'pdf',fn,'start',[0 pi/2 2]);

%% with guessing 
% (doesn't always converge)
% function handle
fn = @(data, mu, sigma, n, g) (1-g).*rd_wrappednormpdf(data, mu, sigma, n) + ...
                                g.*1/(2*pi);

% data to fit
N = 1000;
x = [normrnd(pi/8,pi/4,1,N-N*0.1) -pi + 2*pi.*rand(1,N*0.1)];

% mle
[phat, pci] = mle(x,'pdf',fn,'start',[0 pi/2 2 .05]);

%% with guessing, no bias
% function handle
fn = @(data, sigma, n, g) (1-g).*rd_wrappednormpdf(data, 0, sigma, n) + ...
                                g.*1/(2*pi);

% data to fit
N = 200;
x = [normrnd(0,pi/4,1,N-N*0.1) -pi + 2*pi.*rand(1,N*0.1)];

% mle
[phat, pci] = mle(x,'pdf',fn,'start',[pi/2 2 .05]);

% check fit
[nn xx] = hist(x,30);
y = fn(linspace(-pi,pi,100),phat(1),phat(2),phat(3));
figure
hold on
plot(xx,nn/sum(nn*diff(xx(1:2))))
plot(linspace(-pi,pi,100),y,'r')
    