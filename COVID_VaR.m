clear; close all;
load  Chapter05Data
% whos
%   Date      100x1               801  datetime              
%   Rt        100x1               800  double                
%% Display the Returns
figure(1)
plot(Rt)
xlabel('Time')
ylabel('Returns')
title('Asset Returns')

%% Compute Descriptive Statistics -- See page 116
Stats.Mean=mean(Rt);
Stats.Std=std(Rt);
Stats.Skew=skewness(Rt);
Stats.Kurt=kurtosis(Rt);
Stats.Max=max(Rt);
Stats.Min=min(Rt);
%% Display Histogram with fitted Normal and t-Student -- see page 117
figure(2)
clf
subplot(2,1,1)
histfit(Rt,15)
title('Fit Normal Distribution')
subplot(2,1,2)
histfit(Rt,15, 'tlocationscale' )
title('Fit t-Student Distribution')
%% Display qqplot ( Check functions qqplot, normplot, probplot)
% Is the normality assumption satisfied?
figure()
qqplot(Rt)

%% Test Normality assumption by Jarque-Bera test  - h = jbtest(x)
% Is the normality assumption satisfied?
% You must be able to explain the test to your classmates
h=jbtest(Rt);


%% Probability Levels Corresponding to different values -- See Table 5.2 page 119
mu=Stats.Mean;
sigma=Stats.Std;
normcdf(0.0162,mu,sigma)
alpha=(Rt-mu)./sigma;
[ProbLevel]= normcdf(alpha) ;
[alphaSorted,indexSorted]=sort(alpha);
ProbLevelSorted=ProbLevel(indexSorted);
figure()
clf
plot(alphaSorted, ProbLevelSorted);
xlabel('(Rt-mu)./sigma')
ylabel('Cumulative Probability')
[ProbLevelSorted,alphaSorted]
% test
%% Compute VaR -- Check approximation when hp: mu=0
clear alpha
alpha=0.05;
VaR_1_alpha=mu+norminv(alpha)*sigma
% Variance-Covariance: hp:mu=0 
VaR_1_alpha=norminv(alpha)*sigma


%% Portfolio VaR - Page 140-142
clear
MV=[10,15,20];
Beta=[1.4,1.2 0.8];
Volatility=[0.15,0.12,0.1];
Volatility_Index=0.07;
CorrMat=[1 0.5 0.8;0.5 1 0; 0.8 0 1]; %% CorreMat=ones(3,3);
%CorrMat=eye(3,3)
MV_Portfolio=sum(MV);
Beta_Index=sum(MV.*Beta)/sum(MV)
Virtual=MV.*Beta;
Virtual_Portfolio=MV_Portfolio*Beta_Index

CL=0.99;
Alpha=norminv(CL)

IndividualVaR=Alpha.*Volatility.*MV
PortfolioVaR_IndStocks=sqrt(IndividualVaR*CorrMat*IndividualVaR')
PortfolioVaR_Beta=Virtual_Portfolio*Volatility_Index*Alpha

%% What happens if CorrMat=ones(3,3)?
