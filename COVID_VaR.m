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

%% NOT READY %%%

% mdl1_eu=arima('AR',NaN,'MA',NaN,'Distribution','t','Variance',gjr(1,1));
% WS=22;
% SL=0.05;
% 
% for i=1:length(logReteuro)-WS
%     fit_eu{1,i}=estimate(mdl1_eu,logReteuro(i:i+WS-1));
%     [residuals(:,i),variances(:,i)]=infer(fit_eu{1,i},logReteuro(i:i+WS-1));
%     [muF(i),YMSE(i),sigmaF(i)]=forecast(fit_eu{1,i},1,logReteuro(i:i+WS-1));
%     Parametric_VaR95_eu(i)=-(muF(i)+sigmaF(i)*tinv(SL,fit_eu{1,i}.Distribution.DoF);
% end

%% Parametric VaR across whole sample assuming delta normal

%mu = mean / expected value
mu_eu = mean(logReteuro);
mu_sp = mean(logRetSP500);

%alpha
conf_level=0.95;
alph = norminv(1-conf_level);

%sigma=standard deviation
sigma_eu = std(logReteuro);
sigma_sp = std(logRetSP500);

parametric_var_eu = mu_eu+alph*sigma_eu;
parametric_var_sp = mu_sp+alph*sigma_sp;
%% Parametric VaR with rolling window assuming delta normal

%set alpha
conf_level = 0.95;
alpha = norminv(1-conf_level);

%set window size
WS=22;

%compute stats with rolling window for EU
for i=1:length(logReteuro)-WS
    logRet_eu_insample = logReteuro(i:i+WS-1);
    sigma_eu_in_sample = std(logRet_eu_insample);
    mu_eu_in_sample = mean(logRet_eu_insample);
    VaR_eu_in_sample(i) = alpha*sigma_eu_in_sample+mu_eu_in_sample;
end

%compute stats with rolling window for USA
for c=1:length(logRetSP500)-WS
    logRet_sp_insample = logRetSP500(c:c+WS-1);
    sigma_sp_in_sample = std(logRet_sp_insample);
    mu_sp_in_sample = mean(logRet_sp_insample);
    VaR_sp_in_sample(c) = alpha*sigma_sp_in_sample+mu_sp_in_sample;
end

figure(6)
plot(VaR_eu_in_sample)
hold
bar(logReteuro(1+WS-1:end))
title('Parametric VaR of EUROSTOXX600 using 22 day window')
xlabel('Time')
ylabel('Returns')
legend('VaR','logReturns')

figure(7)
plot(VaR_sp_in_sample)
hold
bar(logRetSP500(1+WS-1:end))
title('Parametric VaR of SP500 using 22 day window')
xlabel('Time')
ylabel('Return')
legend('VaR','logReturns')

%% Parametric VaR with simple moving avg - NOT READY 
%time window for volaitility
%also equals n in simple moving avgs equation
WS_v=7;
%for Paramteric VaR uses same alpha and ws from above

%starting idex needs to be adjusted so we can compute volatitlity in window
%prior
for i=(WS_v+1):length(logReteuro)-WS
    %volatility
    logRet_eu_vol_insample=logReteuro(i-WS_v:i)
    for t=2:length(logRet_eu_vol_insample)
        r_var=logRet_eu_vol_insample(t)-logRet_eu_vol_insample(t-1)
        r_var_sqr(t-1)=r_var*r_var
    end
    %can i do this to add all the elements in r_var_sqr?? What is this data
    %structure called? 
    sigma_vol_eu = sqrt(sum(r_var_sqr)/WS_v-1)
    
%     %parametric VaR
%     %HOW DO WE ADJUST THIS WITH THE VOLATILITY FACTOR?
%     logRet_eu_insample = logReteuro(i:i+WS-1);
%     sigma_eu_in_sample = std(logRet_eu_insample);
%     mu_eu_in_sample = mean(logRet_eu_insample);
%     VaR_eu_in_sample(i) = alpha*sigma_eu_in_sample+mu_eu_in_sample;
end


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

figure()
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

figure(1)
hold on
bar(logReteuro(23:end))
plot(-var_EWMA95_eu, 'r')
hold off

vbt2=varbacktest(logReteuro(23:end),var_EWMA95_eu');
summary(vbt2)
result2 = runtests(vbt2);


%% GARCH codes %%%

model_obj = garch(1,1); % create an obj which define the order of the GARCH model (could be even an ARCH if p=0)
% Estimates model parameters
par_estimates_sp = estimate(model_obj,logRetSP500); 
par_estimates_eu = estimate(model_obj,logReteuro);
% Compute conditional variances time series 
cv_sp = infer(par_estimates_sp,logRetSP500); 
cv_eu = infer(par_estimates_eu,logReteuro);

% Then if we want to compute a forecasting:
k = 10;
sp_cv_forecast = forecast(par_estimates_sp,k,'Y0',logRetSP500);
eu_cv_forecast = forecast(par_estimates_eu,k,'Y0',logReteuro);












