function [VaR] = evt_VaR(n,N_u,u,shape,scale,alpha);
% this is the fuction for the esplicit formula of VaR with extreme value
% theory method (POT)
% the input are:
% n that is a scalar representing the number of the sample
% N_u: representing the number of the peaks over the threshold
% u: the threshold
% shape: the estimated shape parameter of the GPD
% scale: the estimated scale parameter of the GPD
% alpha: confidence level of VaR
VaR=u+(scale/shape)*(((n/N_u*alpha)^(-shape))-1);
end

