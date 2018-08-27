function [ final_price] = sim_trade_pattern_mel_U_ponly(S,tau,theta,sigma,code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A function to simmulate a pattern of trade and then generate a trade
% share matrix and a random sample of final goods prices given Melitz (i.e.
% EKK (2011). 
%
% Note the theta here is the regular one, not AL value 1/theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for goods and countries and the sample size for prices

Ngoods = 50000; % Adjust this number if it is running really slow (not too low though).
Ncntry = length(S);

% Parameters for technologies
eta = sigma; 
markup = eta./(eta-1);

% inveta = 1./(1- eta);
% invNgoods = 1./Ngoods;
% high_price = 1*10^7;

% rand('twister',032878+code)
rng(03281978+code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First thing we want to do is to compute some usefull ojects and then 
% cutoffs. Note I will abstract from the constants here and (hope) that it 
% does not matter.

pi_nn = zeros(Ncntry,1);
cost_cutoff = zeros(Ncntry,1);

for kkk = 1:Ncntry
    
        phi_sum = sum(S.*tau(:,kkk).^(-theta));
        % Compute the stuff in sum across the T_i(w_i*tau_ni)^(-theta)
        
        pi_nn(kkk) = S(kkk)./phi_sum;        
        % Given the S's we can compute the implied pi_nn which is
        % porportional to the ubar_n which relates to the cutoff values
        
        cost_cutoff(kkk) = (pi_nn(kkk)./S(kkk)).^(1/theta);

end

max_goods = max(pi_nn);
% Figure out the country who has the most domestic goods. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw ``u's'' to eventually compute prices in each country.

% m = zeros(Ncntry,Ncntry);
% p_index = zeros(Ncntry,1);

% This may be too complicated, but we are going to have one big price
% array. Each matrix will be indext with the row = goods, coloumn =
% country, then the thrid dimension is a particular importing or consuming
% country. 

% This is just to make clear the divide between a good which is indexed by
% a integer number between 1 and Ngoods and a variety which is indexd by
% the country producer.

price_matrix = zeros(Ngoods, Ncntry);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the new part, 2/17/12. The idea is just have a price matrix, so
% rows are goods, and each coloumn represents the price of the country
% variety for that good. The idea here is to carry this around rather than
% the good by country by variety matrix



for jjj = 1: Ncntry
    
    cgoods = round(Ngoods.*(pi_nn(jjj)./max_goods));
    % Figure out how many domestic enterents there are in a particular
    % country. The max country will have all Ngoods.
    
    s = rand(cgoods,1);
    % Draw uniform variables, these are the little `s's''
    
    u_s = pi_nn(jjj).*s;
    % Now we have to unscale the draws. Again here is where the hope that
    % the constants will not play a role. Check this unscaling if correct.

    % Assighn values in the price matrix, note, I'm going to mix up the
    % values so it is set up where in goods space 1 to Ngoods, some just
    % can be produce domestically, rathern than have 1 to Cgoods produced
    % and all zeros from cgoods tp Ngoods. This may be unesscescary
    
    price_matrix(1:cgoods,jjj) = markup.*(u_s./S(jjj)).^(1/theta);
    
%     price_matrix(:,jjj) = price_matrix(randperm(Ngoods)',jjj);
    
    
    % Multiply by the markup here, then assighn to the price matrix. ;
    % Multiply by the markup here, then assighn to the price matrix. 
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now compute exporting decisions...

final_price = zeros(Ncntry.*Ngoods,Ncntry);
% final_prdcr = zeros(Ncntry.*Ngoods,Ncntry);

for im = 1:Ncntry
    
    carry_prices = zeros(Ngoods,Ncntry);
%     rec_producer = zeros(Ngoods,Ncntry);    
    for ex = 1:Ncntry
        
        if ex == im
            carry_prices(:,im)= price_matrix(:,im);
%             rec_procucer(:,im) = ones(Ngoods,1).*im;
            continue
        end
        
        pth = tau(ex,im).*price_matrix(:,ex) <= markup.*cost_cutoff(im);
        % Find the goods that an exporter can produce lower then the cutoff
        % and enter the market. Note, that 0s in the price matrix will be
        % picked up here, but it will assign the importaing coloumn a zero
        % so it is not big deal to track this.
        
        carry_prices(pth==1,ex) = tau(ex,im).*price_matrix(pth,ex);
%         rec_procucer(:,ex) = ones(Ngoods,1).*ex;
        %No assighn the column for the exporting country, in the importers
        %matrix the price of the good.        
    end           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is now where we will compute the trade flow matrix.

%     flow = carry_prices(:,:)~=0;
    % First figure out the non zero enteries here.    
%     num_parts = zeros(Ncntry,1);
%     
%    for ex = 1:Ncntry
%         
%         num_parts(ex) = sum(carry_prices(flow(:,ex)==1,ex).^(1-eta));
%         % Now for a given importer, for each exporter sum over the prices
%         % to the power (1-eta) which is the numerator component of the CES
%         % expenditure share formula.
%    end
   
%    p_index(im) = (sum(num_parts)).^(-1/(eta-1));
%    % Now the price indes is simply the sum of these components, then taken 
%    % to the power -1/(eta-1),  which I lifted from EKK. Check if this is
%    % correct.
%     
%    m(:,im) = num_parts(:)./(p_index(im)).^(1-eta); 
   % Then the expenditure share is simply equally to the top part computed
   % above divided by the price index taken to the power (1-eta). Again,
   % took it from EKK, seems to match up.
   vvv= carry_prices == 0;
   carry_prices(vvv) = NaN;
   final_price(:,im) = carry_prices(:);
%    final_prdcr(:,im) = rec_procucer(:);
%    final_price(:,im) = ;
end

common_set = sum(isnan(final_price),2)==0; 

% disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now sample of prices. As a first step, the pricing protocoal will
% interpert ICP prices as the geometric average accross country varities
% conditional on a variety existing in all countries. That is we can't have
% a empty entry for a price.



final_price = final_price(common_set,:);
% final_prdcr = final_prdcr(common_set,:);
yyy = 1;




% final_price = nanmean(log(price_matrix(1:20,:,1))')


% rand('twister',02071983+code)
% 
% final_price = zeros(sample,Ncntry);
% 
% for sp = 1:sample
% 
%     keep = randi(Ngoods,1,1);
%     
%     while min(sum(price_matrix(keep,:,:))) == 0
%         keep = randi(Ngoods,1,1);
%     end
%     
%     % These prior llines should insure that only goods where a variety is
%     % consumed by all countries is included.
%     
%     for im = 1:Ncntry
%         p = price_matrix(keep,:,im);        
%         final_price(sp,im) = geomean(p(p~=0));
%         % Now go through each country and compute the non-zero geometric
%         % average across prices.
%     end
%         
% end
% clear price_matrix









