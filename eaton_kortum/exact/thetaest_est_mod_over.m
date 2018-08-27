function [mme] = thetaest_est_mod_over(pricemat,trade_mat,theta,price_error,trade_error,code) %#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the code to compute the estimate of theta
% Size some stuff up
[cntry,ngods]=size(pricemat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add mesurment error.
rand('seed',code)

pricemat = pricemat.*exp(price_error.*randn(cntry,ngods)); 

norm_trdata = trade_mat./(repmat(diag(trade_mat),1,cntry))';

trade_mat = exp(log(norm_trdata) + sqrt(trade_error).*randn(cntry,cntry));
trade_mat(eye(cntry)==1) = 1;
trade_mat = trade_mat./(repmat(sum(trade_mat),cntry,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_p = log(pricemat);
dni = zeros(cntry,cntry);
dni2 = zeros(cntry,cntry);
dni3 = zeros(cntry,cntry);
dni4 = zeros(cntry,cntry);

x5p = round(.65.*ngods + .5);
s5p = round(.75.*ngods + .5);
e0p = round(.80.*ngods + .5);
e5p = round(.85.*ngods + .5);
% n5p = round(.95.*c + .5);


% dni4 = zeros(cntry,cntry);
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
%         num3 = pdiff(h(s5p));
%         num4 = pdiff(h(x5p));
      
        % Compute the mean price difference
        den = mean((pdiff));
        % This is the proxy for the difference in the aggregate price of tradables, i.e. the
        % mean across all prices 
                
        dni(exporter,importer) = num - (den);
        dni2(exporter,importer) = num2 - (den);
%         dni3(exporter,importer) = num3 - (den);
%         dni4(exporter,importer) = num4 - (den);
       
    end
end

% Set up the normalized trade matrix trdx
trdx = trade_mat./repmat(diag(trade_mat),1,length(trade_mat));

% Don't include diagonal entries or (possible) zeros.
diag_trade = eye(cntry);
vvv = (diag_trade(:) == 1);

m1 = log(trdx(~vvv))-((-theta).*(dni(~vvv)));
m2 = log(trdx(~vvv))-((-theta).*(dni2(~vvv)));
% m3 = dni(~vvv).*(log(trdx(~vvv))-((-theta).*(dni(~vvv))));
% m4 = dni2(~vvv).*(log(trdx(~vvv))-((-theta).*(dni2(~vvv))));
% m3 = log(trdx(~vvv))-(-theta).*(dni3(~vvv));
% m4 = log(trdx(~vvv))-(-theta).*(dni4(~vvv));



%Output
mme = [m1,m2];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











