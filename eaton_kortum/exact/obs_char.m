function [dstuff] = obs_char(pricemat,trade_mat,d,border,istraded) %#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes correlations between the price data, trade data,
% and gravity variables.

[r,c]=size(pricemat);
cntry = r;
% istraded = ones(1,c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_p = log(pricemat);
ttt = (istraded == 0);
log_p = log_p(:,~ttt);

dni = zeros(cntry,cntry); dni2 = zeros(cntry,cntry);
tau = zeros(cntry,cntry);
c = 62
e5p = round(.85.*c + .5);
% dni3 = zeros(cntry,cntry); dni4 = zeros(cntry,cntry);

for importer = 1:cntry
    for exporter = 1:cntry
        
        % Compute the price difference
        pdiff = (log_p(importer,:)) - (log_p(exporter,:));
        
        % Sort them, the vector is small, so this is effecient.
        [~, h] = sort(pdiff);
        
        % Now take the max and the 2nd max
        num = pdiff(h(end));
        num2 = pdiff(h(e5p));

        % Compute the mean price difference
        den = mean((pdiff));
        % This is the proxy for the difference in the aggregate price of tradables, i.e. the
        % mean across all prices 
                
        dni(exporter,importer) = num - (den);
        dni2(exporter,importer) = num2 - (den);
        tau(exporter,importer) = num2;

    end
end

% Set up the normalized trade matrix trdx
trdx = trade_mat./repmat(diag(trade_mat),1,length(trade_mat));

% Don't include diagonal entries or (possible) zeros.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup a matrix for country fixed effects to use in the regressions below.

xrow = -ones(1,cntry-1);
eee=eye(cntry-2,cntry-2);

dummy = [-ones(cntry-2,1) eee];
xrow(1,1)=-2;
dummy = [xrow;dummy];

size(dummy);
for i = 1:cntry-2
    xrow = -ones(1,cntry-1);
    dum = [eee(:,1:i), -ones(cntry-2,1), eee(:,i+1:end)];
    xrow(1,i+1) = -2 ;
    dum = [xrow;dum];
    size(dum);
    dummy = [dummy;dum];
end

dummy = [ones(cntry-1,cntry-1) + eye(cntry-1,cntry-1);dummy];

t_trdx = trdx;
ff = (t_trdx(:)~=1);
t_trdx = t_trdx(ff);
d_hat = dni2(ff);
gg = (t_trdx(:)==0 |isnan(d_hat));
dummy = dummy(~gg,:);

vvv = (trdx(:) == 0 | trdx(:) == 1 | isnan(dni2(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then compute correlations...

[cc,pp] = corr((tau(~vvv)),log(d(~vvv)));

[cct,ppt] = corr((tau(~vvv)),log(trdx(~vvv)));

Xmat = [ones(sum(~vvv),1) dummy log(d(~vvv)), border(~vvv)];

[bd,bdint,rsd,~,stats]=regress((dni2(~vvv)),Xmat);

[r c]=size(Xmat); s2 = (1/(r-c)).*sum(rsd.^2);
var_usual = s2*inv(Xmat'*Xmat);
se_usual = diag(var_usual).^.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store and report output

dstuff = [bd(end-1:end)',cc, cct];

if nargout < 1

disp('Correlation and P-value: Distance and Trade Cost')
disp([cc, pp])
disp('Correlation and P-value: Trade and Trade Cost')
disp([cct, ppt])
disp('Distance, Border Elasticity')
disp([bd(end-1:end)])
disp('Distance, Border Confidence Intervals')
disp([bdint(end-1:end,:)])
disp('Standard Error')
disp(se_usual(end-1:end))
disp('R-square from regression')
disp([stats(1)])

end