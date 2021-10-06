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



%% Compute the VaR Using the Historical Simulation Method
% Unlike the normal distribution method, the historical simulation (HS) is
% a nonparametric method. It does not assume a particular 
% distribution of the asset returns. Historical simulation forecasts risk
% by assuming that past profits and losses can be used as the distribution
% of profits and losses for the next period of returns. The VaR "today" is computed as the 
% _p_ th-quantile of the last  _N_  returns prior to "today."

Historical95 = zeros(length(TestWindow),1);
Historical99 = zeros(length(TestWindow),1);

for t = TestWindow
    i = t - TestWindowStart + 1;
    EstimationWindow = t-EstimationWindowSize:t-1;
    X = Returns(EstimationWindow);
    Historical95(i) = -quantile(X,pVaR(1)); 
    Historical99(i) = -quantile(X,pVaR(2)); 
end

figure;
plot(DateReturns(TestWindow),[Historical95 Historical99])
ylabel('VaR')
xlabel('Date')
legend({'95% Confidence Level','99% Confidence Level'},'Location','Best')
title('VaR Estimation Using the Historical Simulation Method')

