function dp_matfile(file, RHO_avg, success_cu_cv, i)

FID = file;
M  = matfile(FID,'Writable',true);
M.RHO_avg(i,1) = RHO_avg;
M.success_cu_cv(i,1) = success_cu_cv;

end