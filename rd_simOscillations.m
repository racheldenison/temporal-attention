% rd_simOscillations.m

t = 0:.01:2;
f = [4:8 10];

soas = [.1 .2 .25 .44 .5];
for i = 1:numel(soas)
    soaIdxs(i) = find(soas(i)==t);
end

a = cos(2*pi*f'*t);

dc = (1:numel(f))'*ones(size(t));

phases = a(:,soaIdxs);


figure
plot(t,a)

figure
plot(t,a+dc)

figure
xoffsets = repmat((0:numel(soas))/500, numel(soas), 1)';
plot((repmat(soas,numel(f),1)+xoffsets)', phases', '.', 'MarkerSize',20)
legend(num2str(f'))
xlabel('soa')
ylabel('phase')
