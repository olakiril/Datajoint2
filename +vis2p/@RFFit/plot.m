function plot(obj,varargin)
params = [];
params = getParams(params,varargin);

a= fetch1(obj,'gauss_fit');

m=a(1:2); 
C=diag(a(3:4)); 
cc=a(5)*sqrt(prod(a(3:4))); 
C(1,2)=cc; 
C(2,1)=cc;

plotGaussRF(m,C,2,params);