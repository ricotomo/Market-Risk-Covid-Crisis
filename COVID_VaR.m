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
%lets watch the distribution -> tails

