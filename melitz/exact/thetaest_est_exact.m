function [mme] = thetaest_est_exact(pricemat,trade_mat,istraded)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the code to compute the estimate of theta EK method for use in
% the SW estimation method.

% Some preliminary stuff.
cntry = length(trade_mat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_p = log(pricemat);
ttt = (istraded == 0);
log_p = log_p(:,~ttt);

[~,c]=size(log_p);

dni = zeros(cntry,cntry);
dni2 = zeros(cntry,cntry);
dni3 = zeros(cntry,cntry);
dni4 = zeros(cntry,cntry);

x5p = round(.65.*c + .5);
s5p = round(.75.*c + .5);
e0p = round(.80.*c + .5);
e5p = round(.85.*c + .5);
n5p = round(.95.*c + .5);


for importer = 1:cntry
    for exporter = 1:cntry
        
        % Compute the price difference
        pdiff = (log_p(importer,:)) - (log_p(exporter,:));
        
        % Sort them, the vector is small, so this is effecient.
        [~, h] = sort(pdiff);

        % Now take the max and the 2nd max
        num = pdiff(h(end));
        num2 = pdiff(h(e5p));
        num3 = pdiff(h(s5p));
%         num4 = pdiff(h(x5p));

      
        % Compute the mean price difference
        den = mean((pdiff));
        % This is the proxy for the difference in the aggregate price of tradables, i.e. the
        % mean across all prices 
                
        dni(exporter,importer) = num - (den);
        dni2(exporter,importer) = num2 - (den);
        dni3(exporter,importer) = num3 - (den);
%         dni4(exporter,importer) = num4 - (den);

    end
end

% Set up the normalized trade matrix trdx
trdx = trade_mat./repmat(diag(trade_mat),1,length(trade_mat));

% Don't include diagonal entries or (possible) zeros.
vvv = (trdx(:) == 1);

mme = [log(trdx(~vvv)), (dni(~vvv)), dni2(~vvv), dni3(~vvv)];

 vvv = isinf(mme);
 mme(vvv) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





