function [mcmcresults, chain, s2chain, sschain] = metropolisHastings(mcmcmodel,data,mcmcparam,mcmcoptions)
%METROPOLISHASTINGS Summary of this function goes here
%   Detailed explanation goes here

%gaussian transition kernel (so it cancels out in the formula as 
% it is symmetric and fwd is equal to backward.



%for gaussian likelihood and gaussian priors. 
%no covariance case -- where is there is no diagonal terms


%prior -- mu -- sigma
for i = 1:length(mcmcparam)
mu_prior(i) = mcmcparam{i,1}{1,5};
sigma_prior_element(i) = mcmcparam{i,1}{1,6};
theta_init(i) =  mcmcparam{i,1}{1,2};
upper_limit(i) = mcmcparam{i,1}{1,4};
lower_limit(i) = mcmcparam{i,1}{1,3};
end
sigma_prior = diag(sigma_prior_element);


theta = theta_init;

%proposal distribution
% no covariance for now

if ~isfield(mcmcoptions,'tuning_factor')
    mcmcoptions.tuning_factor = 0.001; % set to 10 percent;
end

proposal_sigma = sigma_prior*mcmcoptions.tuning_factor ;

acceptance_count = 0; % will keep on updating.

for j = 1:mcmcoptions.nsimu


old_theta = theta;
%generate a proposal theta
%only 1 sample generated.
new_theta =  generate_samples(old_theta,proposal_sigma,upper_limit,lower_limit,1); %generate a new theta from center of old theta
% need to transpose
new_theta = new_theta';

%prior eval
prior_eval_old = gaussian(old_theta,mu_prior,sigma_prior);
prior_eval_new = gaussian(new_theta,mu_prior,sigma_prior);

log_prior_eval_old = log(prior_eval_old);
log_prior_eval_new = log(prior_eval_new);

%LL eval
sigmaLL = mcmcmodel.S20; %keeping it fixed for now


num_pts = size(data,1)*size(data,2);

ss_old = mcmcmodel.ssfun(old_theta,data);
ss_new =  mcmcmodel.ssfun(new_theta,data);

%last two terms will anyway cancel out for fixed theta
% log is ln
LL_old = -ss_old/(2*sigmaLL^2) - num_pts*log(sigmaLL) - (num_pts/2) *log(2*pi);
LL_new = -ss_new/(2*sigmaLL^2) - num_pts*log(sigmaLL) - (num_pts/2) *log(2*pi);

LL_diff(j) = LL_new -LL_old;

log_ratio =  log_prior_eval_new + LL_new - (log_prior_eval_old + LL_old );
log_acceptance_prob =  min(0,log_ratio);
acceptance_prob = exp(log_acceptance_prob);

mcmcresults.prob(j,:) = acceptance_prob;

% a uniform random sample b/w 0 and 1
u = rand;


%accept or reject
if u<= acceptance_prob
    theta =  new_theta;
    acceptance_count = acceptance_count +1;
else
    theta = old_theta;
end

% keep of storing theta
chain(j,:) = theta;

%display 
disp(['Trying sample number  ',num2str(j)])

end





%not implemented yet
mcmcresults.LLdiff = LL_diff;
mcmcresults.acceptance_rate = acceptance_count/mcmcoptions.nsimu;

s2chain = 0;
sschain = 0;


end

