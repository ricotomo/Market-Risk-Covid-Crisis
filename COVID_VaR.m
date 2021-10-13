%% TEAM PROJECT

load('DATASET.mat')

%Calculate Continuous Returns
logRetSP500=tick2ret(pt_SP500,'Method','continuous');
logReteuro=tick2ret(pt_euro,'Method','continuous');

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

% Let's take a look at the distribution -> tails
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
% From the plots we can see negative skewness (bigger losses) and kurtosis (fat tails) -> this confirm the stylized fact of assimetry in gain/losses.
% As we can notice from the graphs above, the t-student distribution fits better for both time series

% Compute moments of distribution for both time series
mean_eu=mean(logReteuro);
mean_SP=mean(logRetSP500);
std_eu=std(logReteuro);
std_SP=std(logRetSP500);
skewness_eu=skewness(logReteuro);
skewness_SP=skewness(logRetSP500);
kurtosis_eu=kurtosis(logReteuro);
kurtosis_SP=kurtosis(logRetSP500);
% Put them together into a table for a better visualization 
table(mean_eu,mean_SP,std_eu,std_SP,skewness_eu,skewness_SP,kurtosis_eu,kurtosis_SP)
% NB: we do expect higher VaR for the SP ptf, since it has more volatility in the considered time horizon than the european one.
% VERIFY THE STYLIZED FACT OF VOLUME/VOLATILITY SIMMETRY  

%% Compute the VaR using the Historical Simulation Method

% Rolling historical VaR of Eurostoxx600
pVaR = [0.05 0.01];
WS=22;

for i=1:length(logReteuro)-WS
    Historical_VaR95_eu(i) = -quantile(logReteuro(i:i+WS-1),pVaR(1)); 
    Historical_VaR99_eu(i) = -quantile(logReteuro(i:i+WS-1),pVaR(2)); 
end

% Rolling historical VaR of S&P500
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
% Both return 1: this indicates that the test rejects the null hypothesis (kurtosis and
% skewness are not 0 -> the normal distribution doesn't represent a good fitting for both time series

%% Compute the VaR using the Parametric Approach

% The parametric VaR assumes that returns and volatility follow a normal
% distribution, but due to the fact that past returns are not a good
% benchmark for the new one, we want to have a really reactive method even
% if we use normal assumption, that's not very realistic so we are going to
% use t-distribution

% Display ACF
figure(4)
subplot(2,2,1)
autocorr(logReteuro)
title('Autocorrelation function logretEU')
subplot(2,2,2)
autocorr(logRetSP500)
title('Autocorrelation function logretSP')
subplot(2,2,3)
autocorr(logReteuro.^2)
title('Autocorrelation function squared logretEU')
subplot(2,2,4)
autocorr(logRetSP500.^2)
title('Autocorrelation function squared logretSP')

% Display PACF
figure(5)
subplot(2,1,1)
parcorr(logReteuro)
title('Partial Autocorrelation function logretEU')
subplot(2,1,2)
parcorr(logRetSP500)
title('Partial Autocorrelation function logretSP')

% As we can see from the graphs, there is autocorrelation among returns,
% and also squared returns, so we have to model firstly the mean taking
% into account some lags dependence and secondly ARCH effects on the
% residuals.
% We can compute the conditional mean by using an ARMA(1,1) model with
% GARCH(1:1) variance of residuals to model autocorr.

%%  Set EWMA to estimate sigma for the Parametric Approach

lambda = 0.95;

% VaR S&P500
histfit(logRetSP500,50,'tlocationscale')
t_score = tinv(0.05,5);
ws = 22;

for i = 2 : (length(logRetSP500) - WS)
    sigma_sp(i) = sqrt((1-lambda) * logRetSP500(i-1)^2 + lambda * var(logRetSP500(i:i+ws-1)));
    var_EWMA95_sp(i) = -t_score*sigma_sp(i);
end

figure(6)
hold on
bar(logRetSP500(23:end))
plot(-var_EWMA95_sp, 'r')
hold off

vbt1=varbacktest(logRetSP500(23:end),var_EWMA95');
summary(vbt)
result1 = runtests(vbt);  

% VaR Eurostoxx600
histfit(logReteuro,50,'tlocationscale')
t_score = tinv(0.05,4);
ws = 22;

for i = 2 : (length(logReteuro) - WS)
    sigma_eu(i) = sqrt((1-lambda) * logReteuro(i-1)^2 + lambda * var(logReteuro(i:i+ws-1)));
    var_EWMA95_eu(i) = -t_score*sigma_eu(i);
end

figure(7)
hold on
bar(logReteuro(23:end))
plot(-var_EWMA95_eu, 'r')
hold off

vbt2=varbacktest(logReteuro(23:end),var_EWMA95_eu');
summary(vbt2)
result2 = runtests(vbt2);

%% EVT not ready
 mdl1_eu=arima('AR',NaN,'Distribution','t','Variance',gjr(1,1));
 WS=500;
 
 for i=1:length(logReteuro)-WS
     fit_eu{1,i}=estimate(mdl1_eu,logReteuro(i:i+WS-1));
     [residuals(:,i),variances(:,i)]=infer(fit_eu{1,i},logReteuro(i:i+WS-1));
     [muF(i),YMSE(i),sigmaF(i)]=forecast(fit_eu{1,i},1,logReteuro(i:i+WS-1));
 end
 tail_fraction=0.05;
 filter_residuals=residuals./sqrt(variances);
 autocorr(filter_residuals(:,1))
 for i=1:length(logReteuro)-WS
     tails{i}=paretotails(filter_residuals(:,i),tail_fraction,1-tail_fraction,'kernel');
 end

[P,Q] = boundary(tails{1});  % cumulative probabilities & quantiles at boundaries
y     = sort(filter_residuals(filter_residuals(:,1) < Q(1),1) - Q(1)); % sort exceedances
plot(y, (cdf(tails{1}, y + Q(1)) - P(1))/P(2))
[F,x] = ecdf(y);                 % empirical CDF
hold on
stairs(x, F, 'r')
grid on












