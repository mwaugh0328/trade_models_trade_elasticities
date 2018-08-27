function [zzz] = monte_carlo_proc(boot)
% This is a program to estimate theta, by Michael Waugh and Ina Simonovska
% for the paper The Elasticity of Trade: Estimates and Evidence. 08/03/12
%
% Note this code computes the estimates when the estimation is
% overidentified using the first and second max

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1 Compute Data moments of interest. 
% Here we will just compute the data moments for you as there are issues
% regarding the sharing of the data set.

load fake_data.mat
[mme]=thetaest_est_exact(pmat_30,tradeshare,istraded);

sample = sum(istraded);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main routine that searches for the theta that mathches the
% observed moments to the moments implied by the model. See est_fun_exact
% for the details.

options = optimset('TolFun',10^-3,'TolX',10^-3,'Display','iter');

mone = mme(:,1);
mtwo = mme(:,2:end);% define parameter first
boot = 092113;

nruns = 12;
nsubs = 200;

tic
[theta_ek, fval_ek] = fminbnd(@(x) est_fun_over(x,mtwo,sample,nruns,nsubs,boot,0),log(3),log(8),options);
toc

disp('Estimate of Theta')
disp(exp(theta_ek))
disp('Test Statistic')
disp(length(mme(:,1)).*fval_ek)

zzz = [exp(theta_ek), length(mme(:,1)).*fval_ek];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
