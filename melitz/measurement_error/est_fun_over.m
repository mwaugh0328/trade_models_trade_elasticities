function loss = est_fun_over(theta,mtwo,sample,Nruns,Nsubs,boot,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This data set has the gravity data and the trade share matrix. 

sig_error = max(theta(2),0);
theta = 1./exp(theta(1));

if flag == 1
    load('../../data/trade_grav_est_30.mat')
else
    load trade_grav_fake.mat
end

% [~, tau, ssd, ~]=gravregasym_sim(tradeshare,distance,b,theta);

[~, tau, ssd, ~]=gravregasym_logd(tradeshare,grav_distance,b,theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [~, tau, ssd, err_var]=gravregimpt_sim(tradeshare,distance,b,theta);
% load fake_data.mat
% [ssd,tau] = tau_est_exact(pmat_30,tradeshare,1./theta,istraded);

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
[ndata, nmoments] = size(mtwo);
nmoments = 3;
Ncntry = length(ssd);
Nobvs = Ncntry.^2 - Ncntry;

d_hat = log(grav_distance(eye(Ncntry)~=1));

boot_L = 617429;

record = zeros(Nobvs,Nruns,nmoments);
second_part = zeros(nmoments,nmoments,Nruns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now the simmulatio is performed. Included is the
% .mex and the normal .m file. The .mex runs much faster. 
%
% Given the trade flows, mhat and prices final_price, it passes this to the
% gen moments routine. See this for more details.
%
% Finally, I use the parfor command, this distributes each computation to a
% different core. This obviously speeds things up a lot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if Nruns ~= 1
    parfor runs = 1:Nruns

    [final_price] = sim_trade_pattern_mel_U_ponly(exp(ssd),tau,1./theta,sigma,runs+boot);

    
    [record(:,runs,:)]  = gen_moments(final_price,d_hat,sig_error,Nsubs,sample,nmoments,runs+boot);

    end

else

        runs = 1;
        
        [final_price] = sim_trade_pattern_mel_U_ponly(exp(ssd),tau,1./theta,sigma,runs+boot);
        
        [record(:,runs,:)]  = gen_moments(final_price,d_hat,sig_error,Nsubs,sample,nmoments,runs+boot);

end

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Record stuff and then compute the difference between the observed moments
% % and the simmulated moments.
% 
s_moments = reshape(nanmean(record,2),Nobvs,nmoments);



m1 = mtwo(:,1);
% m2 = (mtwo(:,1)-mean(m1)).^2;
m3 = mtwo(:,2);
% m4 = mtwo(:,3);
% % m4 = (mtwo(:,2)-mean(m3)).^2;
m5 = (mtwo(:,1)-mean(m1)).*(d_hat-mean(d_hat)); % This computes the covariance.
% m6 = (mtwo(:,2)-mean(m3)).*(d_hat-mean(d_hat));
% m3 = mone - ((-1./theta).*mtwo(:,3));
% m4 = mone - ((-1./theta).*mtwo(:,4));

d_moments = [m1,m3,m5];

y_theta = mean(d_moments-s_moments)';

% second_part = (1./Nruns).*sum(second_part,3);
W = eye(nmoments);
% loss = y_theta'*y_theta;

% W = cov(d_moments-s_moments);
loss = y_theta'*(W^-1)*y_theta;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % if Nruns ~= 1
%     
% else
%     loss = y_theta'*y_theta;
% end


