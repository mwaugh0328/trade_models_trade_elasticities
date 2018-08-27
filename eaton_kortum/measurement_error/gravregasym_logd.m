function [tssdmat, rtausd, ssd, err_var]=gravregasym_logd(trdshrs,distance,border,theta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program runs the inital gravity regression that yeilds the S's and
% the taus for a given theta.  
%
% Note to Yasutora, this is the code you need to modify to introduce a free
% trade dummy varible. The two key things are (1) put it in the gravity
% regression below and (2) reconstruct the implied trade costs given the
% estimatiaes. I will flag these belwo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This just constructs a bunch of matricies to have country fixed effects
% and exporter fixed effects in Waugh (2010). Note I have code if you want
% to have importer fixed effects as in EK(2002)

N = length(trdshrs);

xrow = -ones(1,N-1);
eee=eye(N-2,N-2);

dummy = [-ones(N-2,1) eee];
xrow(1,1)=-2;
dummy = [xrow;dummy];

size(dummy);
for i = 1:N-2
    xrow = -ones(1,N-1);
    dum = [eee(:,1:i), -ones(N-2,1), eee(:,i+1:end)];
    xrow(1,i+1) = -2 ;
    dum = [xrow;dum];
    size(dum);
    dummy = [dummy;dum];
end

dummy = [ones(N-1,N-1) + eye(N-1,N-1);dummy];

asym = eye(N-1,N-1);
for i = 1:N-1
    ass = eye(N-1,N-1);
    ass(i,:) = [];
    ass = [-ones(1,N-1);ass];
    asym = [asym;ass];
end

dest(1:N-1,1:N-1)=-1;
for i = 1:N-1
    int = zeros(N-1,N-1);
    int(:,i) = ones(N-1,1);
    dest = [dest;int];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are zeros in the trade data, so our approach is to take them out,
% then when we construct the implied trade costs below, we use the
% estimates to compute the trade costs in the entries where there are
% zeros. This just does some adjustments here.
%
% Note regarding zeros, in the liturature they seem to make a big deal out
% of this. My experince it is not, simple OLS is perfectly fine.

trdata = trdshrs;
trdata = trdata./(repmat(diag(trdata),1,N))';
hh = (trdata(:)~=0 & trdata(:)~=1);
qq = (trdata(:)~=0);
ff = (trdata(:)~=1);
trdata = trdata(ff);
gg = (trdata(:)~=0);
trdata = trdata(gg);
dummy = dummy(gg,:);
asym = asym(gg,:);
dest = dest(gg,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mdist = .6213711*distance;
mdist_nz = mdist(ff);
mdist = mdist(hh);
for i = 1:length(distance(hh))
    ddd  = mdist(i);
    dmat(i,:) = double([(ddd < 375) (375<=ddd & ddd<750) (750<=ddd& ddd<1500) (1500<=ddd& ddd<3000)...
            (3000<=ddd& ddd<6000) (6000<=ddd)]);
end

for i = 1:length(distance(ff))
    ddd  = mdist_nz(i);
    dmat_nz(i,:) = double([(ddd < 375) (375<=ddd & ddd<750) (750<=ddd& ddd<1500) (1500<=ddd& ddd<3000)...
            (3000<=ddd& ddd<6000) (6000<=ddd)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is the regression
border_nz = border(ff);
border = border(hh);

% Here you want to put in the free trade dummy.
[bsd,~,~,~,statssd] = regress(log(trdata),[dummy asym  log(mdist) ones(length(mdist),1) border]);

err_var = statssd(end);

ssd = [-sum(bsd(1:N-1));bsd(1:N-1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct Fitted Values (again, make sure to include new free trade dummy)
ysd =exp([dummy asym log(mdist) ones(length(mdist),1) border ]*bsd);

msdmat = zeros(N^2,1);
msdmat(hh) = ysd;
msdmat(~ff) = 1;
msdmat = reshape(msdmat,N,N);
sdc = sum(msdmat);
sdc = repmat(sdc,N,1);
tssdmat = msdmat./sdc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Trade Costs (again, make sure to adjust for the free trade dummy)


bsasym = bsd(N:end-3);
bsasym = [-sum(bsasym);bsd(N:end-3)];
beasym = exp(-theta.*bsasym)-1;

befsd = exp(-theta.*bsd(end,:))-1;

beasym=(beasym(:,:)*ones(1,N));
beasym = beasym(ff);

dist_elas = bsd(end-2,:);

rtausd = zeros(N^2,1);

rtausd(~ff) = 1;
rtausd(ff) = exp(-theta.*(bsd(end-1,:) + log(mdist_nz).*dist_elas)).*(1+border_nz.*befsd).*(1+beasym);
rtausd = reshape(rtausd,N,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



