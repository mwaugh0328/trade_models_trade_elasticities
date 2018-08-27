function loss = gen_fake_date(theta,sample,boot)

sig_error = 0;
theta = 1./theta(1);

load estimation_mat_30.mat
[~, tau, ssd, err_var]=gravregasym_logd(tradeshare,distance,b,theta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifies values for the simmulation. 
% Here sigma plays no real role.
%
% Nruns is the numer of times a new trade matrix is simmulated with the
% associated prices. The trade flows vary little across the simmulations,
% hence the small number. 
%
% For each trade flow price data set, prices are sampled. This is done 100 
% times per simmulation of the trade flows. The number of prices sampled is
% 50, the same as in the EK data set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 1.5;
 
[mhat, final_price] = sim_trade_pattern_ek_mex(exp(ssd),tau,1./theta,sigma,boot);
    

[Ngoods,~] = size(final_price);
rng(12071940+boot)

keep = randi(Ngoods,sample,1);

mhat = trade_add_error(mhat, err_var,12071940+boot);

final_price_tilde = (final_price(keep,:))';
pmat_30 = final_price_tilde.*lognrnd(0,sig_error,length(mhat),length(final_price_tilde));  
tradeshare = mhat;

istraded = ones(1,length(pmat_30));

save fake_data pmat_30 tradeshare istraded

save trade_grav_fake tradeshare b distance 

loss = 1;
        
        
        
        