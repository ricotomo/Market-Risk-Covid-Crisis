function [VaR] = evt_VaR(n,N_u,u,shape,scale,alpha);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
VaR=u+(scale/shape)*(((n/N_u*alpha)^(-shape))-1);
end

