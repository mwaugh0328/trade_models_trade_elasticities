clear
load ek_results

Nsims = 100;
results = zeros(Nsims,3);

for xxx = 1:Nsims
    
    gen_fake_date(theta,sample,05012014+xxx);
    
    results(xxx,:)=monte_carlo_proc(0328);
    
    save standard_errors_new xxx results
    
end
