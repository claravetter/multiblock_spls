function cv_mbspls_outer_mean_deflate(analysisfolder)

resultfile = [analysisfolder, 'final_results/preliminary_result_mean_boot.mat'];
load(resultfile)

%% deflation
matrices_deflated = input.Xs; 
for num_m=1:size(input.Xs,2)
    weights{num_m} = output.final_parameters{1,3}; % weight vector for matrix num_m
    matrices_deflated{num_m} = cv_proj_def_single(input.Xs{num_m}, weights{1}{num_m});
end
input.Xs = matrices_deflated; 
%%
input.covariates{1,4} = [];
input.covariates{1,6} = [];
input.covariates_names{1,4} = [];
input.covariates_names{1,6} = [];
%%
datafile = [analysisfolder, 'final_results/preliminary_result_mean.mat'];
input.datafile = datafile; 
input.CV = output.CV; 
input.save_CV = 0;
save(datafile, 'input', 'setup', 'output' )

end