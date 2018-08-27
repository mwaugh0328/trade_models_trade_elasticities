function moments = gen_moments(rec_low_price,d_hat,sig_error,Nsubs,sample,nmoments,boot)
% This code computes the moments from simmulated trade flow and price data

code = 072279+boot;
[Ngoods,Ncntry] = size(rec_low_price);
rng(code,'twister')
cn = Ncntry.^2 - Ncntry;

moments = zeros(cn,Nsubs,nmoments);

for sub_runs = 1:Nsubs
        
        % Add disturbances to the trade shares, the err_var is the st.dev
        % to the error term in the gravity regression
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Draw the number of goods that we want
            
        keep = randi(Ngoods,sample,1);

        final_price_tilde = (rec_low_price(keep,:))';
        
%         m = trade_add_error(mhat, err_var,007+sub_runs+boot);
    
        moments(:,sub_runs,:) = thetaest_est_mod_D(final_price_tilde, d_hat,sig_error,Nsubs+code);
        
        
%         moments(:,sub_runs,:) = thetaest_est_mod_over(final_price_tilde, m, theta);

        %         moments(sub_runs,:) = thetaest_est_mod_med_mean_mex(final_price_tilde, m);
%         moments(sub_runs,:) = thetaest_est_mod_over_mex(final_price_tilde, m);
        % using the .mex provides alot of speed up. I will include the
        % regular file as well.
    
end

moments = mean(moments,2);
vvv = isinf(moments);
moments(vvv) = NaN;

