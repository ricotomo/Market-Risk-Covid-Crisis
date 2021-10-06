<<<<<<< Updated upstream
%% TEAM PROJECT

%Calculate Continuos Returns
logRetSP500=log(pt_SP500(2:end)./pt_SP500(1:end-1))
logReteuro=log(pt_euro(2:end)./pt_euro(1:end-1))

%% Display the Returns
figure(1)
plot(logRetSP500)
xlabel('Time')
ylabel('Returns')
title('Returns SP500')

figure(2)
plot(logReteuro)
xlabel('Time')
ylabel('Returns')
title('Returns STOXX600')
%lets watch the distribution -> tails
=======
%% TEAM PROJECT 
% 1) Historical simulation 
>>>>>>> Stashed changes
