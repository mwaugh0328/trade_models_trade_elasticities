function [mme] = thetaest_est_mod_D(pricemat,distance,price_error,code) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the code to compute the estimate of theta
% Size some stuff up
[cntry,ngods]=size(pricemat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add mesurment error.
%randn('seed',code);
rng(code)

pricemat = pricemat.*exp(price_error.*randn(cntry,ngods)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_p = log(pricemat);
dni = zeros(cntry,cntry);
dni2 = zeros(cntry,cntry);
dni3 = zeros(cntry,cntry);

e5p = round(.85.*ngods + .5);
s5p = round(.75.*ngods + .5);

% Now we are going to go through by country pair, compute the price
% differences, take max, second max, etc.
for importer = 1:cntry
    for exporter = 1:cntry
        
        % Compute the price difference
        pdiff = (log_p(importer,:)) - (log_p(exporter,:));
        
        % Sort them, the vector is small, so this is effecient.
        [g, h] = sort(pdiff);

        % Now take the max and the 2nd max
        num = pdiff(h(end));
        num2 = pdiff(h(e5p));
        num3 = pdiff(h(s5p));

      
        % Compute the mean price difference
        den = mean((pdiff));
        % This is the proxy for the difference in the aggregate price of tradables, i.e. the
        % mean across all prices 
                
        dni(exporter,importer) = num - (den);
        dni2(exporter,importer) = num2 - (den);
        dni3(exporter,importer) = num3 - (den);

       
    end
end

% Set up the normalized trade matrix trdx
% trdx = trade_mat./repmat(diag(trade_mat),1,length(trade_mat));

% Don't include diagonal entries or (possible) zeros.
diag_trade = eye(cntry);
vvv = (diag_trade(:) == 1);

m1 = dni(~vvv);
% m2 = (dni(~vvv)-mean(m1)).^2;
m3 = dni2(~vvv);
m4 = dni3(~vvv);
% % m4 = (dni2(~vvv)-mean(m3)).^2;
%m5 = ((dni(~vvv)-mean(m1)).*(distance-mean(distance)))./(std(m1).*std(distance));
% m6 = (dni2(~vvv)-mean(m3)).*(distance-mean(distance));

%Output
mme = [m1];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











