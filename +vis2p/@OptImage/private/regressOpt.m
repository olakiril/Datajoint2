function [B, R2, Fp, dof] = regressOpt(X, G, addDoF)
% Linear regression
%
% INPUTS:
%   X      TxN matrix of signals, where T is the # samples and N is # of signals
%   G      TxK matrix of predictors where K is the number of modeled predictors
%   addDoF additional degrees of freedom that are not included in the
%          design matrix, e.g. when multiple design matrices 
%
% OUTPUTS:
%   B    regression cofficients of the model predictors.
%   R2   real-valued relative response amplitude, e.g. dF/F
%   Fp   p-value of each predictor computed from t-distribution
%   dof  degrees of freedom in residual based on autocorrelation

assert(size(G,1)==size(X,1))
assert(isreal(X))

B = (G/(G'*G))'*X;              % regression coefficients
R2 = sum(X.^2);                 % initial variance
xf = abs(fft(X));               % power spectrum
dof = sum(xf).^2./sum(xf.^2);   % degrees of freedom in original signal
X = X - G*B;                    % residual
R2 = 1-sum(X.^2)./R2;           % R-squared
dof1 = size(B,1)+addDoF;        % degrees of freedom in model
Fp = 1-fcdf(R2.*dof/dof1, dof1, dof);   % p-value of the F distribution
end