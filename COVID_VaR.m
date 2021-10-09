%% TEAM PROJECT
load('DATASET.mat')
%Calculate Continuous Returns
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
% Compute Statistics
mean_eu=mean(logReteuro);
mean_SP=mean(logRetSP500);
std_eu=std(logReteuro);
std_SP=std(logRetSP500);
skewness_eu=skewness(logReteuro);
skewness_SP=skewness(logRetSP500);
kurtosis_eu=kurtosis(logReteuro);
kurtosis_SP=kurtosis(logRetSP500);
% Main statistics of the two indexes
table(mean_eu,mean_SP,std_eu,std_SP,skewness_eu,skewness_SP,kurtosis_eu,kurtosis_SP)
% NB: we do expect higher VaR for the SP ptf because it has more
% volatility(in this time horizon) than the european one.
%% Compute the VaR Using the Historical Simulation Method
% Rolling historical VaR of Eurostoxx
pVaR = [0.05 0.01];
WS=22;
for i=1:length(logReteuro)-WS
    Historical_VaR95_eu(i) = -quantile(logReteuro(i:i+WS-1),pVaR(1)); 
    Historical_VaR99_eu(i) = -quantile(logReteuro(i:i+WS-1),pVaR(2)); 
end
% Rolling historical VaR of SP500

for i=1:length(logRetSP500)-WS
    Historical_VaR95_SP(i) = -quantile(logRetSP500(i:i+WS-1),pVaR(1)); 
    Historical_VaR99_SP(i) = -quantile(logRetSP500(i:i+WS-1),pVaR(2)); 
end

figure(3)
plot(Dates_eu(24:end),Historical_VaR95_eu)
hold
plot(Dates_SP(24:end),Historical_VaR95_SP,'r')
xlabel('Time')
ylabel('VaR at 95%')
title('Historical VaR at 95% US vs EU')
legend('EU','US')
%% Testing the Normality assumption by Jarque-Bera test
h_eu=jbtest(logReteuro);
h_sp=jbtest(logRetSP500);
%both return 1
%this indicates that the test rejects the null hypothesis (kurtosis and
%skewness are not 0 i.e. the distribution does not have a good normal fit)

%% Compute the VaR Using the parametric approach
%the parametric VaR assumes that returns and volatility follow a normal
%distribution

%mu = mean / expected value
mu_eu=mean(logReteuro)
mu_sp=mean(logRetSP500)

%alpha
conf_level=0.95;
alpha=norminv(1-conf_level)

%sigma=standard deviation
sigma_eu=std(logReteuro)
sigma_sp=std(logRetSP500)

parametric_var_eu=mu_eu+alpha*sigma_eu



%% Compute the VaR Using the Parametric Method









