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
hold on
plot(Dates_SP(24:end),Historical_VaR95_SP,'r')
xlabel('Time','FontSize', 13)
ylabel('VaR at 95%','FontSize', 13)
title('Historical VaR at 95% US vs EU','FontSize', 15)
legend('EU VaR','US VaR', 'FontSize', 13)
hold off

figure(4)
subplot(1,2,1)
hold on 
bar(Dates_eu(502:end),logReteuro(501:end), 'r')
plot(Dates_eu(502:end),-Historical_VaR95_eu(479:end),'b-')
xlabel('Time','FontSize', 15)
ylabel('VaR at 95%','FontSize', 15)
title('Historical VaR w/ returns EU','FontSize',16)
legend('Returns','Historical VaR 95%','FontSize', 13)
hold off
subplot(1,2,2)
hold on
bar(Dates_SP(502:end),logRetSP500(501:end), 'r')
plot(Dates_SP(502:end),-Historical_VaR95_SP(479:end),'g-')
xlabel('Time', 'FontSize',15)
ylabel('VaR at 95%','FontSize', 15)
title('Historical VaR w/ returns US','FontSize',16)
legend('Returns','Historical VaR 95%','FontSize', 13)
hold off

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

% Eventually we decide to choose dof of freedom equal to 3, since lower than 3 the variance wouldn't be finite 

% VaR S&P500

% 5% Confidence level 
t_score = tinv(0.05,3);
ws = 22;

for i = 2 : (length(logRetSP500) - ws)
    sigma_sp(i) = sqrt((1-lambda) * logRetSP500(i-1)^2 + lambda * var(logRetSP500(i:i+ws-1)));
    var_EWMA95_sp(i) = -t_score*sigma_sp(i);
end

% 10% Confidence level 
t_score = tinv(0.1,3);
ws = 22;

for i = 2 : (length(logRetSP500) - ws)
    sigma_sp(i) = sqrt((1-lambda) * logRetSP500(i-1)^2 + lambda * var(logRetSP500(i:i+ws-1)));
    var_EWMA90_sp(i) = -t_score*sigma_sp(i);
end

% VaR Eurostoxx600

% 5% Confidence level
t_score = tinv(0.05,3);
ws = 22;

for i = 2 : (length(logReteuro) - ws)
    sigma_eu(i) = sqrt((1-lambda) * logReteuro(i-1)^2 + lambda * var(logReteuro(i:i+ws-1)));
    var_EWMA95_eu(i) = -t_score*sigma_eu(i);
end

% 10% Confidence level
t_score = tinv(0.1,3);
ws = 22;

for i = 2 : (length(logReteuro) - ws)
    sigma_eu(i) = sqrt((1-lambda) * logReteuro(i-1)^2 + lambda * var(logReteuro(i:i+ws-1)));
    var_EWMA90_eu(i) = -t_score*sigma_eu(i);
end

figure(5)
plot(Dates_eu(24:end),var_EWMA95_eu)
hold on
plot(Dates_SP(24:end),var_EWMA95_sp,'r')
xlabel('Time','FontSize', 15)
ylabel('VaR at 95%','FontSize', 15)
title('Parametric VaR at 95% EU vs US','FontSize', 15)
legend('EU VaR','US VaR', 'FontSize', 15)
hold off

figure(6)
subplot(1,2,1)
hold on 
bar(Dates_eu(502:end),logReteuro(501:end), 'r')
plot(Dates_eu(502:end),-var_EWMA95_eu(479:end),'b-')
xlabel('Time','FontSize', 15)
ylabel('VaR at 95%','FontSize', 15)
title('Parametric VaR w/ returns EU','FontSize',16)
legend('Returns','EWMA VaR 95%','FontSize', 13)
hold off
subplot(1,2,2)
hold on
bar(Dates_SP(502:end),logRetSP500(501:end), 'r')
plot(Dates_SP(502:end),-var_EWMA95_sp(479:end),'b-')
xlabel('Time', 'FontSize',15)
ylabel('VaR at 95%','FontSize', 15)
title('Parametric VaR w/ returns US','FontSize',16)
legend('Returns','EWMA VaR 95%','FontSize', 13)
hold off

%% Display ACF
figure(7)
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

% EVT STOXX, WS 500 
% Estimate AR(1) and EGARCH(1;1) and forecasts for each WS of sigma and mu
 mdl1_eu=arima('AR',NaN,'Distribution','t','Variance',gjr(1,1));
 WS=500;
 
 for i=1:length(logReteuro)-WS
     fit_eu{1,i}=estimate(mdl1_eu,logReteuro(i:i+WS-1));
     [residuals_eu(:,i),variances_eu(:,i)]=infer(fit_eu{1,i},logReteuro(i:i+WS-1));
     [muF_eu(i),YMSE_eu(i),sigmaF_eu(i)]=forecast(fit_eu{1,i},1,logReteuro(i:i+WS-1));
 end
 
% VaR STOXX 90%
% Fitting filter residuals tail with GDP 
 tail_fraction_eu=0.1;
 filter_residuals_eu=residuals_eu./sqrt(variances_eu);
 
 figure(8)
 autocorr(filter_residuals_eu(:,1))
 
 for i=1:length(logReteuro)-WS
     tails_eu_90{i}=paretotails(filter_residuals_eu(:,i),tail_fraction_eu,1-tail_fraction_eu,'kernel');
 end
 
 for i=1:length(logReteuro)-WS
    [P(:,i),Q(:,i)]=boundary(tails_eu_90{i});
    params(i,:)=lowerparams(tails_eu_90{i});
 end
 
% Computing Dynamic VaR 90% 
for i=1:length(logReteuro)-WS
    N_u(i)=length(filter_residuals_eu(filter_residuals_eu(:,i) < Q(1,i),1));
[quantile_90_evt_eu(i)]=evt_VaR(WS,N_u(i),Q(1,i),params(i,1),params(i,2),0.1);
 VaR_90_evt_eu(i)=muF_eu(1,i)+sqrt(sigmaF_eu(1,i))*(-quantile_90_evt_eu(i));
end

% Fitted tail distribution vs empirical (for the first Window Size)
x_sorted_eu=sort(filter_residuals_eu(:,1));
x_eu=x_sorted_eu(1:N_u(1));
F_y_eu=cdf(tails_eu_90{1},x_eu);
F_x_eu=F_y_eu/0.1;
[F_em_eu,x_em_eu] = ecdf(x_eu);   % empirical CDF

figure(9)
plot(x_eu,F_x_eu); % estimated
hold on
stairs(x_em_eu, F_em_eu, 'r.','MarkerSize',20)
grid on
xlabel('Residuals','FontSize',13)
ylabel('Cumulative Probability','FontSize',13)
title('Estimated vs empirical tail STOXX','FontSize',16)

% VaR STOXX 95%
% Fitting filter residuals tail with GDP 
 tail_fraction_eu=0.05;
 
 for i=1:length(logReteuro)-WS
     tails_eu_95{i}=paretotails(filter_residuals_eu(:,i),tail_fraction_eu,1-tail_fraction_eu,'kernel');
 end
 
 for i=1:length(logReteuro)-WS
    [P(:,i),Q(:,i)]=boundary(tails_eu_95{i});
    params(i,:)=lowerparams(tails_eu_95{i});
 end

for i=1:length(logReteuro)-WS
    N_u(i)=length(filter_residuals_eu(filter_residuals_eu(:,i) < Q(1,i),1));
[quantile_95_evt_eu(i)]=evt_VaR(WS,N_u(i),Q(1,i),params(i,1),params(i,2),0.05);
 VaR_95_evt_eu(i)=muF_eu(1,i)+sqrt(sigmaF_eu(1,i))*(-quantile_95_evt_eu(i));
end

% EVT S&P, WS 500  
% Estimate AR(1) and EGARCH(1;1) and forecasts for each WS of sigma and mu
 mdl1_sp=arima('AR',NaN,'Distribution','t','Variance',gjr(1,1));
 WS=500;
 
 for i=1:length(logRetSP500)-WS
     fit_sp{1,i}=estimate(mdl1_sp,logRetSP500(i:i+WS-1));
     [residuals_sp(:,i),variances_sp(:,i)]=infer(fit_sp{1,i},logRetSP500(i:i+WS-1));
     [muF_sp(i),YMSE_sp(i),sigmaF_sp(i)]=forecast(fit_sp{1,i},1,logRetSP500(i:i+WS-1));
 end
 
 % VaR S&P500 90%
 % Fitting filter residuals tail with GDP 
 tail_fraction_sp=0.1;
 filter_residuals_sp=residuals_sp./sqrt(variances_sp);
 
 figure(10)
 autocorr(filter_residuals_sp(:,1))
 
 for i=1:length(logRetSP500)-WS
     tails_sp_90{i}=paretotails(filter_residuals_sp(:,i),tail_fraction_sp,1-tail_fraction_sp,'kernel');
 end
 
  for i=1:length(logRetSP500)-WS
    [P(:,i),Q(:,i)]=boundary(tails_sp_90{i});
    params(i,:)=lowerparams(tails_sp_90{i});
  end
  
 % Computing dynamic EVT VaR 90%  
for i=1:length(logRetSP500)-WS
    N_u(i)=length(filter_residuals_sp(filter_residuals_sp(:,i) < Q(1,i),1));
[quantile_90_evt_sp(i)]=evt_VaR(WS,N_u(i),Q(1,i),params(i,1),params(i,2),0.1);
 VaR_90_evt_sp(i)=muF_sp(1,i)+sqrt(sigmaF_sp(1,i))*(-quantile_90_evt_sp(i));
end

% Fitted tail distribution vs empirical (for the first Window Size)
x_sorted_sp=sort(filter_residuals_sp(:,1));
x_sp=x_sorted_sp(1:N_u(1));
F_y_sp=cdf(tails_sp_90{1},x_sp);
F_x_sp=F_y_sp/0.1;
[F_em_sp,x_em_sp] = ecdf(x_sp);   % empirical CDF

figure(11)
plot(x_sp,F_x_sp); % estimated
hold on
stairs(x_em_sp, F_em_sp, 'r.','MarkerSize',20)
grid on
xlabel('Residuals','FontSize',13)
ylabel('Cumulative Probability','FontSize',13)
title('Estimated vs empirical tail S&P500','FontSize',16)

% VaR S&P500 95% 
 tail_fraction_sp=0.05;
 
  for i=1:length(logRetSP500)-WS
     tails_sp_95{i}=paretotails(filter_residuals_sp(:,i),tail_fraction_sp,1-tail_fraction_sp,'kernel');
  end
 
 for i=1:length(logRetSP500)-WS
    [P(:,i),Q(:,i)]=boundary(tails_sp_95{i});
    params(i,:)=lowerparams(tails_sp_95{i});
 end
 
for i=1:length(logRetSP500)-WS
    N_u(i)=length(filter_residuals_sp(filter_residuals_sp(:,i) < Q(1,i),1));
[quantile_95_evt_sp(i)]=evt_VaR(WS,N_u(i),Q(1,i),params(i,1),params(i,2),0.05);
 VaR_95_evt_sp(i)=muF_sp(1,i)+sqrt(sigmaF_sp(1,i))*(-quantile_95_evt_sp(i));
end

%% Display the results with returns for both US and EU at 90% and 95%

figure(12)
subplot(2,2,1)
plot(Dates_eu(502:end),-VaR_90_evt_eu)
hold
bar(Dates_eu(502:end),logReteuro(501:end),'r')
xlabel('Time')
ylabel('VaR EVT 90% and  log-returns EU')
title('VaR EVT 90% vs log-returns EU')
subplot(2,2,2)
plot(Dates_SP(502:end),-VaR_90_evt_sp)
hold
bar(Dates_SP(502:end),logRetSP500(501:end),'r')
xlabel('Time')
ylabel('VaR EVT 90% and  log-returns SP')
title('VaR EVT 90% vs log-returns SP')
subplot(2,2,3)
plot(Dates_eu(502:end),-VaR_95_evt_eu)
hold
bar(Dates_eu(502:end),logReteuro(501:end),'r')
xlabel('Time','FontSize',13)
ylabel('VaR EVT 95% and  log-returns EU','FontSize',13)
title('VaR EVT 95% vs log-returns EU','FontSize',16)
subplot(2,2,4)
plot(Dates_SP(502:end),-VaR_95_evt_sp)
hold
bar(Dates_SP(502:end),logRetSP500(501:end),'r')
xlabel('Time','FontSize',13)
ylabel('VaR EVT 95% and  log-returns SP','FontSize',13)
title('VaR EVT 95% vs log-returns SP','FontSize',16)

%% ANSWER TO THE FIRST RESEARCH QUESTION (MARKET HYPOTHESIS)

figure(13)
subplot(1,2,1)
plot(Dates_eu(502:end),VaR_90_evt_eu)
hold
plot(Dates_SP(502:end),VaR_90_evt_sp)
xlabel('Time')
ylabel('VaR 90%')
title('VaR 90% STOXX vs SP500')
legend('STOXX','S&P500')
subplot(1,2,2)
plot(Dates_eu(502:end),VaR_95_evt_eu)
hold
plot(Dates_SP(502:end),VaR_95_evt_sp)
xlabel('Time')
ylabel('VaR 95%')
title('VaR 95% STOXX vs SP500')
legend('STOXX','S&P500')

%% ANSWER TO THE SECOND RESEARCH QUESTION (RISK HYPOTHESIS)

% We can see that during the covid crisis the evt is the best model with 90%
figure(14)
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

% With 95%
figure(15)
subplot(1,2,1)
plot(Dates_eu(502:end-40),-VaR_95_evt_eu(1:end-40),'g-')
hold
bar(Dates_eu(502:end-40),logReteuro(501:end-40))
plot(Dates_eu(502:end-40),-Historical_VaR95_eu(479:end-40),'b-')
plot(Dates_eu(502:end-40),-var_EWMA95_eu(479:end-40),'m')
xlabel('Time')
ylabel('VaR at 95%')
title('Comparison of models during Covid crisis STOXX')
legend('VaR 95% EVT','Returns','VaR 95% Historical','VaR 95% EWMA')
subplot(1,2,2)
plot(Dates_SP(502:end-40),-VaR_95_evt_sp(1:end-40),'g-')
hold
bar(Dates_SP(502:end-40),logRetSP500(501:end-40))
plot(Dates_SP(502:end-40),-Historical_VaR95_SP(479:end-40),'b-')
plot(Dates_SP(502:end-40),-var_EWMA95_sp(479:end-40),'m')
xlabel('Time')
ylabel('VaR at 95%')
title('Comparison of models during Covid crisis S&P500')
legend('VaR 95% EVT','Returns','VaR 95% Historical','VaR 95% EWMA')

%% BACKTESTING

% At 90% confidence level 

% Historical back test US 
vbt_his_sp=varbacktest(logRetSP500(23:end),Historical_VaR90_SP','VaRLevel',.90);
summary(vbt_his_sp)
result_his_sp = runtests(vbt_his_sp); 
% Historical back test EU 
vbt_his_eu=varbacktest(logReteuro(23:end),Historical_VaR90_eu','VaRLevel',.90);
summary(vbt_his_eu)
result_his_eu = runtests(vbt_his_eu);

% Parametric back test US
vbt_par_sp=varbacktest(logRetSP500(23:end),var_EWMA90_sp','VaRLevel',.90);
summary(vbt_par_sp)
result_par_sp = runtests(vbt_par_sp);  
% Parametric back test EU
vbt_par_eu=varbacktest(logReteuro(23:end),var_EWMA90_eu','VaRLevel',.90);
summary(vbt_par_eu)
result_par_eu = runtests(vbt_par_eu);

% EVT back test SP
vbt_evt_sp=varbacktest(logRetSP500(501:end),VaR_90_evt_sp(1:end)','VaRLevel',.90);
summary(vbt_evt_sp)
result_evt_sp = runtests(vbt_evt_sp);
% EVT back test EU
vbt_evt_eu=varbacktest(logReteuro(501:end),VaR_90_evt_eu(1:end)','VaRLevel',.90);
summary(vbt_evt_eu)
result_evt_eu = runtests(vbt_evt_eu);

% At 95% confidence level 

% Historical back test US 
vbt_his_sp_95=varbacktest(logRetSP500(23:end),Historical_VaR95_SP');
summary(vbt_his_sp_95)
result_his_sp_95 = runtests(vbt_his_sp_95); 
% Historical back test EU 
vbt_his_eu_95=varbacktest(logReteuro(23:end),Historical_VaR95_eu');
summary(vbt_his_eu_95)
result_his_eu_95 = runtests(vbt_his_eu_95);

% Parametric back test US
vbt_par_sp_95=varbacktest(logRetSP500(23:end),var_EWMA95_sp');
summary(vbt_par_sp_95)
result_par_sp_95 = runtests(vbt_par_sp_95);  
% Parametric back test EU
vbt_par_eu_95=varbacktest(logReteuro(23:end),var_EWMA95_eu');
summary(vbt_par_eu_95)
result_par_eu_95 = runtests(vbt_par_eu_95);

% EVT back test SP
vbt_evt_sp_95=varbacktest(logRetSP500(501:end),VaR_95_evt_sp(1:end)');
summary(vbt_evt_sp_95)
result_evt_sp_95 = runtests(vbt_evt_sp_95);
% EVT back test EU
vbt_evt_eu_95=varbacktest(logReteuro(501:end),VaR_95_evt_eu(1:end)');
summary(vbt_evt_eu_95)
result_evt_eu_95 = runtests(vbt_evt_eu_95);

figure(16)
plot(Dates_eu(1:500),Q(1,1).*ones(1,500),'r')
hold
plot(Dates_eu(1:500),filter_residuals_eu(:,1)','b.','MarkerSize',10)
title('Peaks over the thresholds','FontSize',16)
xlabel('Time','FontSize',13)
ylabel('Filter Residuals','FontSize',13)







