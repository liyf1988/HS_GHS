%This file gives a working example of the HSGHS function

clear all
n = 100; p = 50; q = 25;

%Random generator for coefficients and precision matrix
num_nonzero = p*q/5;
B = zeros(1,p*q);
ind = randsample(p*q, num_nonzero);
B(ind) = 0.5+rand([1 num_nonzero])*1.5;
B(ind(1:num_nonzero/2)) = -B(ind(1:num_nonzero/2));
B = reshape(B,p,q);
%True precision matrix has an AR1 structure
Omega_true = toeplitz([1,0.45,zeros(1,q-2)]);
Sigma_true = eye(q)/Omega_true;

%Generate design matrix (n by p)
X = mvnrnd(zeros(n,p), eye(p));

%Generate observed values
E = mvnrnd(zeros(n,q), Sigma_true);
Y = X*B + E;

%HS-GHS estimate
burnin = 1000; nmc = 5000;
[beta_save,lambda_sq_save,tau_sq_save,...
omega_save,lambda_G_sq_save,tau_G_sq_save] = HSGHS(X,Y,burnin,nmc,eye(q));

%%%Get posterior mean
beta_mean = mean(beta_save,3);
omega_mean = mean(omega_save,3);
