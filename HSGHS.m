% Yunfan Li, August 2018
% Sampling for coefficients and precision in multivariate regression, using HS and GHS

function [beta_save,lambda_sq_save,tau_sq_save,...
omega_save,lambda_G_sq_save,tau_G_sq_save] = HSGHS(X,Y,burnin,nmc,Sigma_init)
% Input:
%     X: n*p matrix of observed predictors
%     Y: n*q matrix of observed responses
%     burnin, nmc: number of MCMC burnins and saved samples
%     Sigma_init: initial value of covariance matrix; need to be symmetric and pos-definite q*q matrix

% Output:
%     beta_save: p by q by nmc matrices of saved posterior samples of
%     coefficients
%     lambda_sq_save: p*q by nmc vector of saved samples of lambda
%     squared (local tuning parameter)
%     tau_sq_save: 1 by nmc vector of saved samples of tau squared (global
%     tuning parameter)
%     omega_save: q by q by nmc matrices of saved posterior samples of
%     precision matrix
%     lambda_sq_G_save: q*(q-1)/2 by nmc vector of saved samples of lambda
%     squared (local tuning parameter), for precision matrix
%     tau_sq_G_save: 1 by nmc vector of saved samples of tau squared (global
%     tuning parameter), for precision matrix

[n,q] = size(Y);
p = size(X,2);

beta_save = zeros(p,q,nmc);
lambda_sq_save = zeros(p,q,nmc);
tau_sq_save = zeros(nmc,1);
omega_save = zeros(q,q,nmc);
lambda_G_sq_save = zeros(q*(q-1)/2,nmc);
tau_G_sq_save = zeros(nmc,1);

ind_all = zeros(q-1,q);
for i = 1:q
       if i==1  
       ind = [2:q]'; 
      elseif i==p
       ind = [1:q-1]'; 
      else
       ind = [1:i-1,i+1:q]';
       end
       
       ind_all(:,i) = ind;
end

% set initial values
beta = zeros(p,q); lambda_sq = ones(p*q,1); nu = ones(p*q,1);
tau_sq = 1; xi = 1; 
Omega = eye(q); Sigma = Sigma_init;
Lambda_G_sq(1:q,1:q) = 1; Nu_G(1:q,1:q) = 1; tau_G_sq = 1; xi_G = 1;

for iter = 1:(burnin+nmc)
          
    if(mod(iter,1000)==0)
		fprintf('iter = %d \n',iter);
    end

%%% reshape and scale data
	Omega_chol = chol(Omega);
	tmp = Omega_chol*Y';
	y_scaled = tmp(:);
	X_scaled = kron(X,Omega_chol);
	
%%% sample beta
	lambda_tau = lambda_sq.*tau_sq;
	u = normrnd(0,1,[p*q 1]).*sqrt(lambda_tau);
	delta = normrnd(0,1,[n*q 1]);
	v = X_scaled*u + delta;
	w = (X_scaled*diag(lambda_tau)*(X_scaled')+eye(n*q))\(y_scaled-v);
	beta = u + diag(lambda_tau)*(X_scaled'*w);
%%% sample lambda_sq and nu
    rate = beta.^2/(2*tau_sq)+1./nu;
    lambda_sq = 1./gamrnd(1,1./rate);    % random inv gamma with shape=1, rate=rate
    nu = 1./gamrnd(1,1./(1+1./lambda_sq));    % random inv gamma with shape=1, rate=1+1/lambda_sq_12
%%% sample tau_sq and xi	
    rate = 1/xi + sum(beta.^2./(2*lambda_sq));
    tau_sq = 1/gamrnd((p*q+1)/2, 1/rate);    % inv gamma w/ shape=(p*(p-1)/2+1)/2, rate=rate
    xi = 1/gamrnd(1,1/(1+1/tau_sq));    % inv gamma w/ shape=1, rate=1+1/tau_sq

%%% subtract mean, for selected subsample
	B_sample = reshape(beta,[q p]);
	B_sample = B_sample';
	Y_fitted = X*B_sample;
	Y_res = Y - Y_fitted;
	S = Y_res'*Y_res;

%%% sample Sigma and Omega=inv(Sigma)
    for i = 1:q
      ind = ind_all(:,i);     
      Sigma_11 = Sigma(ind,ind); sigma_12 = Sigma(ind,i);
      sigma_22 = Sigma(i,i);
      s_21 = S(ind,i); s_22 = S(i,i);
      lambda_G_sq_12 = Lambda_G_sq(ind,i); nu_G_12 = Nu_G(ind,i);
      %% sample gamma and beta_G
      gamma = gamrnd((n/2+1),2/s_22);    % random gamma with shape=n/2+1, rate=s_22/2
      inv_Omega_11 = Sigma_11 - sigma_12*sigma_12'/sigma_22;
      inv_C = s_22*inv_Omega_11+diag(1./(lambda_G_sq_12*tau_G_sq));
      inv_C_chol = chol(inv_C);
      mu_i = -inv_C\s_21;
      beta_G = mu_i+ inv_C_chol\randn(q-1,1);
      omega_12 = beta_G; omega_22 = gamma + beta_G'*inv_Omega_11*beta_G;
      %% sample lambda_G_sq and nu_G
      rate = omega_12.^2/(2*tau_G_sq)+1./nu_G_12;
      lambda_G_sq_12 = 1./gamrnd(1,1./rate);    % random inv gamma with shape=1, rate=rate
      nu_G_12 = 1./gamrnd(1,1./(1+1./lambda_G_sq_12));    % random inv gamma with shape=1, rate=1+1/lambda_sq_12
      %% update Omega, Sigma, Lambda_G_sq, Nu_G
      Omega(i,ind) = omega_12; Omega(ind,i) = omega_12;
      Omega(i,i) = omega_22;
      temp = inv_Omega_11*beta_G;
      Sigma_11 = inv_Omega_11 + temp*temp'/gamma;
      sigma_12 = -temp/gamma; sigma_22 = 1/gamma;
      Sigma(ind,ind) = Sigma_11; Sigma(i,i) = sigma_22;
      Sigma(i,ind) = sigma_12; Sigma(ind,i) = sigma_12;
      Lambda_G_sq(i,ind) = lambda_G_sq_12;
	  Lambda_G_sq(ind,i) = lambda_G_sq_12;
      Nu_G(i,ind) = nu_G_12; Nu_G(ind,i) = nu_G_12;
    end
    
%%% sample tau_sq and xi
    omega_vector = Omega(tril(true(size(Omega)),-1));
    lambda_G_sq_vector = Lambda_G_sq(tril(true(size(Lambda_G_sq)),-1));
    rate = 1/xi_G + sum(omega_vector.^2./(2*lambda_G_sq_vector));
    tau_G_sq = 1/gamrnd((q*(q-1)/2+1)/2, 1/rate);    % inv gamma w/ shape=(p*(p-1)/2+1)/2, rate=rate
    xi_G = 1/gamrnd(1,1/(1+1/tau_G_sq));    % inv gamma w/ shape=1, rate=1+1/tau_sq

%%% save B, Omega, shrinkage parameters
    if iter > burnin           
		beta_save(:,:,iter-burnin) = B_sample;
                tmp = reshape(lambda_sq,[q p]); tmp = tmp';
		lambda_sq_save(:,:,iter-burnin) = tmp;
		tau_sq_save(iter-burnin) = tau_sq;
		omega_save(:,:,iter-burnin) = Omega;
        lambda_G_sq_save(:,iter-burnin) = lambda_G_sq_vector;
        tau_G_sq_save(iter-burnin) = tau_G_sq;
	end

end

end