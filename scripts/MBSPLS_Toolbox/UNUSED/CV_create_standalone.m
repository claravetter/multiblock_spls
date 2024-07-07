%  script CV_create_standalone

 
spls_standalone_version = 'Dev_2024';


mkdir('/opt/SPLS/')
cd('/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox')
eval(['mcc -m  /volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox_standalone/' spls_standalone_version '...
          '-d /volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox_standalone/main']);
        
test = system('qsub /volume/DP_FEF/ScrFun/ScriptsRepository/DP_SPLS_standalone.sh');

