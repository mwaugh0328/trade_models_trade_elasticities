function yyyy = trade_add_error(p_trade, err_var,code)

N = length(p_trade);

p_trdata = p_trade./(repmat(diag(p_trade),1,N))';
ff = (p_trdata(:)~=1);

% randn('seed',12131943+code)
rng(12131943+code)

msdmat = zeros(N^2,1);
msdmat(ff) = exp(log(p_trdata(ff)) + sqrt(err_var).*randn(length(p_trdata(ff)),1));
msdmat(~ff) = 1;
msdmat = reshape(msdmat,N,N);
sdc = sum(msdmat);
sdc = repmat(sdc,N,1);
yyyy = msdmat./sdc;