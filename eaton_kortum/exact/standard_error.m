clear
load ek_results

Nsims = 500;
results = zeros(Nsims,2);



for xxx = 1:Nsims
    
    gen_fake_date(theta,sample,2565587+xxx);
    
    results(xxx,:)=monte_carlo_proc(0328);
    
    save standard_errors xxx results
    
end

