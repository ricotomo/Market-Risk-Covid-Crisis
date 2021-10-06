%% TEAM PROJECT
load('DATASET.mat')
%Calculate Continuos Returns
logRetSP500=log(pt_SP500(2:end)./pt_SP500(1:end-1));
logReteuro=log(pt_euro(2:end)./pt_euro(1:end-1));

%% Display the Returns
figure(1)
subplot(2,1,1)
plot(Dates_SP(2:end),logRetSP500)
xlabel('Time')
ylabel('LogReturns')
title('LogReturns SP500')
subplot(2,1,2)
plot(Dates_eu(2:end),logReteuro,'r')
xlabel('Time')
ylabel('LogReturns')
title('LogReturns STOXX600')
% We can see volatility clustering 
%lets watch the distribution -> tails
figure(2)
subplot(2,2,1)
histfit(logReteuro)
xlabel('logRet')
ylabel('Frequency')
title('Fitting Normal Distribution STOXX600')
subplot(2,2,2)
histfit(logRetSP500)
xlabel('logRet')
ylabel('Frequency')
title('Fitting Normal Distribution SP500')
subplot(2,2,3)
histfit(logReteuro,100,'tlocationscale')
xlabel('logRet')
ylabel('Frequency')
title('Fitting t-Distribution STOXX600')
subplot(2,2,4)
histfit(logRetSP500,100,'tlocationscale')
xlabel('logRet')
ylabel('Frequency')
title('Fitting t-Distribution SP500')
% we can see negative skewness (bigger losses) and kurtosis (fat tails),
% t-student fits better 





