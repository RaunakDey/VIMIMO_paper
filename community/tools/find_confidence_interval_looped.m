function [S_min,S_max,V_min,V_max,S_median,V_median] = find_confidence_interval_looped(chain,transient_id_new,mcmcpars,confidence_limit,model, pars_without_nan)
%FIND_CONFIDENCE_INTERVALS : by Raunak Dey: gives confidence limits, max + min; inputs:
%chain, %of intervals, 
%   Detailed explanation goes here

chain_to_be_used = chain(transient_id_new:end,:);
[length_chain_to_be_used,~]=size(chain_to_be_used);



tvec = 0:0.05:15.75; % for better viz

i=1;
for j = 1:10:length_chain_to_be_used
    
    %pars_from_dist = @(chain_to_be_used) chain_to_be_used;
    pars_used = update_pars(pars2,chain_to_be_used(i,:),mcmcpars);
    [~,S_recon(:,:,i),V_recon(:,:,i),D_recon(:,:,i)] = simulate_ode(model,pars_used,tvec,pars_used.S0,pars_used.V0); % mcmc parameter set
    i=i+1
end



for i = 1:model.NH
     [~,~,l] = size(S_recon)
    S_reconstructed_i = reshape(S_recon(:,i,:),length(tvec),l);
    S_min(:,i) = quantile(S_reconstructed_i',(1-confidence_limit)/2);
    S_max(:,i) = quantile(S_reconstructed_i',1-(1-confidence_limit)/2);
    S_median(:,i) = median(S_reconstructed_i');
    clear S_reconstructed_i;
end



for i=1:model.NV

    V_reconstructed_i = reshape(V_recon(:,i,:),length(tvec),l);
    V_min(:,i) = quantile(V_reconstructed_i',(1-confidence_limit)/2);
    V_max(:,i) = quantile(V_reconstructed_i',1-(1-confidence_limit)/2);
    V_median(:,i) = median(V_reconstructed_i');
    clear V_reconstructed_i;

end


end





