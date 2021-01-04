%Run simulations for Figure 3
%Requires results from runSimTable1.m and analyseSimTable1.m

%%Bootstrap the median MSE
simoptions=3*simoptions+1;

t=tempname();
mkdir(t);
disp(['temporary space for parcluster: ', t]);

% some cleanup
c=parcluster('local');
c.JobStorageLocation=t;
delete(c.Jobs);

WavePath
addpath(strcat(pwd,'/FISTA'))
addpath(strcat(pwd,'/code'))
addpath(strcat(pwd,'/simfunc'))

options=astro_setparams();

simoptions=simoptions-1;
i=mod(simoptions,24)+1;
simoptions=floor(simoptions/24);
nM=simoptions+1;

M=2^(nM+7)

load(strcat('data_and_results/EOsim',num2str(M),'.mat'))
load(strcat('data_and_results/rsummary_powerlaw_',num2str(M),'.mat'))
	
clear result
curresult=rsummary(i);

if(i<=6) 
    type='blockssym'; %cosmoblocks
elseif(i<=12) 
    type='f4sym'; %cosmo2
elseif(i<=18) 
    type='truef'; %cosmo1
elseif(i<=24) 
    type='f4'; %cosmo2asymmetric
end

if(M==256) 
    nboot=48
    parpool(12)
elseif(M==512)
    nboot=24
    parpool(4)
end

options.ps=curresult.ps;
options.nitr=500;
options.showplot=0;
options.parallel=1;
%options.parallel=0;
options.center=0;

[sol,F]= astrobootWS(curresult.F,M/2,E,0,nboot,options);

save(strcat('data_and_results/rsummary_diag_powerlaw_',num2str(M),type,'ps',num2str(options.ps),'type',num2str(curresult.type*2),'nboot',num2str(nboot),'.mat'),'curresult','sol');

