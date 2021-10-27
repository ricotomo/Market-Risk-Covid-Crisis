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
pVaR = [0.05 0.1];
WS=22;

for i=1:length(logReteuro)-WS
    Historical_VaR95_eu(i) = -quantile(logReteuro(i:i+WS-1),pVaR(1)); 
    Historical_VaR90_eu(i) = -quantile(logReteuro(i:i+WS-1),pVaR(2)); 
end

% Rolling historical VaR of S&P500
for i=1:length(logRetSP500)-WS
    Historical_VaR95_SP(i) = -quantile(logRetSP500(i:i+WS-1),pVaR(1)); 
    Historical_VaR90_SP(i) = -quantile(logRetSP500(i:i+WS-1),pVaR(2)); 
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

%%  Set EWMA to estimate sigma for the Parametric Approach

lambda = 0.95;

% Find the dof that fits best the data (fitting a t-student distribution)
% SP timeseries
pd_sp = fitdist(logRetSP500, 'tlocationscale');
disp(pd_sp.nu)

% euro timeseries
pd_euro = fitdist(logReteuro, 'tlocationscale');
disp(pd_euro.nu)

% Eventually we decide to choose dof of freedom equal to 3/4, since lower than 3 the variance wouldn't be finite 

% VaR S&P500
% 5%
histfit(logRetSP500,50,'tlocationscale')
t_score = tinv(0.05,3);
ws = 22;

for i = 2 : (length(logRetSP500) - ws)
    sigma_sp(i) = sqrt((1-lambda) * logRetSP500(i-1)^2 + lambda * var(logRetSP500(i:i+ws-1)));
    var_EWMA95_sp(i) = -t_score*sigma_sp(i);
end
% 10%
t_score = tinv(0.1,3);
ws = 22;

for i = 2 : (length(logRetSP500) - ws)
    sigma_sp(i) = sqrt((1-lambda) * logRetSP500(i-1)^2 + lambda * var(logRetSP500(i:i+ws-1)));
    var_EWMA90_sp(i) = -t_score*sigma_sp(i);
end
figure(6)
hold on
bar(logRetSP500(23:end))
plot(-var_EWMA95_sp, 'r')
hold off

% VaR Eurostoxx600
% 5%
histfit(logReteuro,50,'tlocationscale')
t_score = tinv(0.05,3);
ws = 22;

for i = 2 : (length(logReteuro) - ws)
    sigma_eu(i) = sqrt((1-lambda) * logReteuro(i-1)^2 + lambda * var(logReteuro(i:i+ws-1)));
    var_EWMA95_eu(i) = -t_score*sigma_eu(i);
end

% 10%
histfit(logReteuro,50,'tlocationscale')
t_score = tinv(0.1,3);
ws = 22;

for i = 2 : (length(logReteuro) - ws)
    sigma_eu(i) = sqrt((1-lambda) * logReteuro(i-1)^2 + lambda * var(logReteuro(i:i+ws-1)));
    var_EWMA90_eu(i) = -t_score*sigma_eu(i);
end
figure(7)
hold on
bar(logReteuro(23:end))
plot(-var_EWMA95_eu, 'r')
hold off

%% Display ACF
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

% As we can see from the graphs, there is autocorrelation among returns,
% and also squared returns, so we have to model firstly the mean taking
% into account some lags dependence and secondly ARCH effects on the
% residuals.
% We can compute the conditional mean by using an AR(1) model with
% GARCH(1:1) variance of residuals to model autocorr.
%% EXTREME VALUE THEORY
% EVT 500 WS STOXX
% estimate AR(1) and EGARCH(1;1) and forecasts for each WS of sigma and mu
 mdl1_eu=arima('AR',NaN,'Distribution','t','Variance',gjr(1,1));
 WS=500;
 
 for i=1:length(logReteuro)-WS
     fit_eu{1,i}=estimate(mdl1_eu,logReteuro(i:i+WS-1));
     [residuals_eu(:,i),variances_eu(:,i)]=infer(fit_eu{1,i},logReteuro(i:i+WS-1));
     [muF_eu(i),YMSE_eu(i),sigmaF_eu(i)]=forecast(fit_eu{1,i},1,logReteuro(i:i+WS-1));
 end
 % fitting filter residuals tail with GDP 
 tail_fraction_eu=0.1;
 filter_residuals_eu=residuals_eu./sqrt(variances_eu);
 figure(6)
 autocorr(filter_residuals_eu(:,1))
 for i=1:length(logReteuro)-WS
     tails_eu{i}=paretotails(filter_residuals_eu(:,i),tail_fraction_eu,1-tail_fraction_eu,'kernel');
 end
 
% Adjusting residuals (exceedances) and the associated  cumulative probability
for i=1:length(logReteuro)-WS
    x_sorted_eu(:,i)=sort(filter_residuals_eu(:,i));
 
end
  x_eu=x_sorted_eu(1:50,:);
for i=1:length(logReteuro)-WS  
    [F_y_eu(:,i)]=cdf(tails_eu{i},x_eu(:,i));
    
end
F_x_eu=F_y_eu/0.1;

% Fitted tail distribution vs empirical (for the first Window Size)
[F_em_eu,x_em_eu] = ecdf(x_eu(:,1));   % empirical CDF
figure(7)
plot(x_eu(:,1),F_x_eu(:,1)); % estimated
hold on
stairs(x_em_eu(:,1), F_em_eu(:,1), 'r.')
grid on
xlabel('Residuals')
ylabel('Probability/Frequency')
title('Estimated vs empirical tail STOXX')

% computing VaR with dynamic EVT     
 for i=1:length(logReteuro)-WS
    [~,idx5_eu(i)] = min(abs(F_x_eu(:,i)-0.1));
    quantile_90_evt_eu(i)=x_eu(idx5_eu(i),i);
    VaR_90_evt_eu(i)=muF_eu(1,i)+sqrt(sigmaF_eu(1,i))*(-quantile_90_evt_eu(i));
 end

% EVT 500 WS S&P
% estimate AR(1) and EGARCH(1;1) and forecasts for each WS of sigma and mu
 mdl1_sp=arima('AR',NaN,'Distribution','t','Variance',gjr(1,1));
 WS=500;
 
 for i=1:length(logRetSP500)-WS
     fit_sp{1,i}=estimate(mdl1_sp,logRetSP500(i:i+WS-1));
     [residuals_sp(:,i),variances_sp(:,i)]=infer(fit_sp{1,i},logRetSP500(i:i+WS-1));
     [muF_sp(i),YMSE_sp(i),sigmaF_sp(i)]=forecast(fit_sp{1,i},1,logRetSP500(i:i+WS-1));
 end
 % fitting filter residuals tail with GDP 
 tail_fraction_sp=0.1;
 filter_residuals_sp=residuals_sp./sqrt(variances_sp);
 figure(8)
 autocorr(filter_residuals_sp(:,1))
 for i=1:length(logRetSP500)-WS
     tails_sp{i}=paretotails(filter_residuals_sp(:,i),tail_fraction_sp,1-tail_fraction_sp,'kernel');
 end
 
 % Adjusting residuals (exceedances) and the associated  cumulative probability
for i=1:length(logRetSP500)-WS
    x_sorted_sp(:,i)=sort(filter_residuals_sp(:,i));
 
end
  x_sp=x_sorted_sp(1:50,:);
for i=1:length(logRetSP500)-WS  
    [F_y_sp(:,i)]=cdf(tails_sp{i},x_sp(:,i));
    
end
F_x_sp=F_y_sp/0.1;

% Fitted tail distribution vs empirical (for the first Window Size)
[F_em_sp,x_em_sp] = ecdf(x_sp(:,1));   % empirical CDF
figure(9)
plot(x_sp(:,1),F_x_sp(:,1)); % estimated
hold on
stairs(x_em_sp(:,1), F_em_sp(:,1), 'r.')
grid on
xlabel('Residuals')
ylabel('Probability/Frequency')
title('Estimated vs empirical tail SP500')

% computing VaR with dynamic EVT     
 for i=1:length(logRetSP500)-WS
    [~,idx5_sp(i)] = min(abs(F_x_sp(:,i)-0.1));
    quantile_90_evt_sp(i)=x_sp(idx5_sp(i),i);
    VaR_90_evt_sp(i)=muF_sp(1,i)+sqrt(sigmaF_sp(1,i))*(-quantile_90_evt_sp(i));
 end
%% Display the results with returns for both
figure(10)
subplot(1,2,1)
plot(Dates_eu(502:end),-VaR_90_evt_eu)
hold
bar(Dates_eu(502:end),logReteuro(501:end))
xlabel('Time')
ylabel('VaR EVT and  log-returns EU')
title('VaR EVT vs log-returns EU')
subplot(1,2,2)
plot(Dates_SP(502:end),-VaR_90_evt_sp)
hold
bar(Dates_SP(502:end),logRetSP500(501:end))
xlabel('Time')
ylabel('VaR EVT and  log-returns SP')
title('VaR EVT vs log-returns SP')
%% ANSWER TO THE FIRST RESEARCH QUESTION (MARKET HYPOTHESIS)
figure(11)
plot(Dates_eu(502:end),VaR_90_evt_eu)
hold
plot(Dates_SP(502:end),VaR_90_evt_sp)
xlabel('Time')
ylabel('VaR 90%')
title('VaR 95% STOXX vs SP500')
legend('STOXX','S&P500')
%% ANSWER TO THE SECOND RESEARCH QUESTION (RISK HYPOTHESIS)
% we can see that during the covid crisis the evt is the best model
figure(12)
subplot(1,2,1)
plot(Dates_eu(502:end),-VaR_90_evt_eu,'g-')
hold
bar(Dates_eu(502:end),logReteuro(501:end))
plot(Dates_eu(502:end),-Historical_VaR90_eu(479:end),'b-')
plot(Dates_eu(502:end),-var_EWMA90_eu(479:end),'m')
xlabel('Time')
ylabel('VaR at 90%')
title('Comparison of models during Covid crisis STOXX')
legend('VaR 90% EVT','Returns','VaR 90% Historical','VaR 90% EWMA')
subplot(1,2,2)
plot(Dates_SP(502:end),-VaR_90_evt_sp,'g-')
hold
bar(Dates_SP(502:end),logRetSP500(501:end))
plot(Dates_SP(502:end),-Historical_VaR90_SP(479:end),'b-')
plot(Dates_SP(502:end),-var_EWMA90_sp(479:end),'m')
xlabel('Time')
ylabel('VaR at 90%')
title('Comparison of models during Covid crisis S&P500')
legend('VaR 90% EVT','Returns','VaR 90% Historical','VaR 90% EWMA')


%% BACKTESTING
% Historical back test US 
vbt1=varbacktest(logRetSP500(23:end),Historical_VaR95_SP');
summary(vbt1)
result1 = runtests(vbt1); 
% Historical back test EU 
vbt2=varbacktest(logReteuro(23:end),Historical_VaR95_eu');
summary(vbt2)
result2 = runtests(vbt2);

% Parametric back test US
vbt3=varbacktest(logRetSP500(23:end),var_EWMA95_sp');
summary(vbt3)
result3 = runtests(vbt3);  
% Parametric back test EU
vbt4=varbacktest(logReteuro(23:end),var_EWMA95_eu');
summary(vbt4)
result4 = runtests(vbt4);

% EVT back test SP
vbt_evt_sp=varbacktest(logRetSP500(501:end),VaR_90_evt_sp','VaRLevel',.90);
summary(vbt_evt_sp)
% EVT back test EU
vbt_evt_eu=varbacktest(logReteuro(501:end),VaR_90_evt_eu','VaRLevel',.90);
summary(vbt_evt_eu)











