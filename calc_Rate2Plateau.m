function [r,p,m,v] = calc_Rate2Plateau(data,limit,n,vars)
% r = rate, as in Nstrides to plateau
% p = upper and lower limits of plateau value
% m = mean of plateau value
% v = variability of plateau data - limit +/- mean

% unpacking vars
plt = vars(1);
group = vars(2);
condition = num2str(vars(3));
subject = num2str(vars(4));
day = num2str(vars(5));

gttl = [{'80%'},{'100%'},{'120%'},{'FPa'},{'STa'},{'S2a'}];
% data is a vector of one individuals parameter, limit is the value the
% parameter must reach in #SDs from the final n strides

% excluding data points greater than 1.5 SDs from the model fit
g = fittype('a+b*exp(-c*x)');
rMax = 0.1;
ldata = length(data);
xdata = (1:ldata)';
startPoints = [[ones(size(xdata)), -exp(-xdata)]\data; 1];
fit1 = fit(xdata,data,g,'Start',startPoints,...
    'lower',[-Inf -Inf 0],'upper',[Inf Inf rMax],'Robust','on');
fdata = feval(fit1,xdata);

% limiting data to 300 strides, esp for day 2
len = ldata;
if ldata > 300
    len = 300;
end
Res = fdata(1:len) - data(1:len);
I = abs(Res) > 1.5*std(data(1:len));
dataOut = data(~I); % data without outliers

% mean of last n strides
m = nanmean(dataOut(end-n:end));

% variability of last n strides
% v = limit*nanstd(data(end-n:end));
% variability of residuals
v = limit*nanstd(Res(round((1/3)*len):end));

% upper and lower limits
ul = m+v; ll = m-v;

% plateau limits
p = [{'upperlimit'},{'lowerlimit'};
    {ul},{ll}];

% strides within the upper and lower limits specified by limit*sd
inLimit = data < ul & data > ll;

% summing inLimit to find bins of 30 strides with 30 values within upper
% and lower limits
for ii = 1:length(data)-n
    s = sum(inLimit(ii:ii+(n-1)));
    if s == 30
        break
    end
end

r = ii;

if plt == 1
    figure
    plot(fit1)
    hold on
    plot(data,'.k')
    plot([r r],[-1 1],'k')
    plot([0 600],[ll ll],'r')
    plot([0 600],[ul ul],'r')
    title([gttl{group} ' subject: ' subject ' adapt: ' condition ' day: ' day])
    ylabel('L_{sym})')
    xlabel('strides')
    xlim([0 xdata(end)])
    text(10,0.8,['stride_{plateau}: ' num2str(r)],'FontSize',12,'FontWeight','Bold')
    legend off
end
    
    




