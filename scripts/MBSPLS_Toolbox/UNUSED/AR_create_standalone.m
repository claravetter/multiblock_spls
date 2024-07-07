%  script AR_create_standalone
% Fun = function to create a standalone from
% Dir = directory of the subfunctions

mkdir('/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox_standalone')
cd('/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox')
eval(['mcc -m  /volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox_standalone/dp_spls_standalone.m '...
          '-d /volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox_standalone/main']);
        
test = system('qsub /volume/DP_FEF/ScrFun/ScriptsRepository/DP_SPLS_standalone.sh');

